## ==========================================================================================
## inferPath Algorithm: Ligand-Receptor-Target-Regulon Signal Pathway Inference
## ==========================================================================================
## Method Overview:
## To infer molecular pathways of signal regulation between cells/domains, inferPath algorithm:
## 1. Ligands from source domains bind to receptors on target domains
## 2. Signal transduction through intracellular intermediate nodes activates transcriptional regulators
## 3. Regulators initiate expression of downstream target genes
## 4. Activity of transcription factors and target genes is denoted as regulons
## 
## Workflow:
## Step 1: CellChat calculates enriched ligand-receptor (L-R) pairs between source and target domains
## Step 2: NicheNet predicts ligand activities and infers ligand-target regulatory links
## Step 3: SCENIC regulons are integrated to connect targets to transcription factors
## Step 4: Construct ligand-receptor-target-regulon signaling network
## ==========================================================================================


args <- commandArgs(trailingOnly = T)

scale.ratio.sig <- as.numeric(args[1]) #0.8
topdegnum <- as.integer(args[2])  # 1000
stages <- args[3]


library(CellChat)
library(patchwork)
library(dplyr)
library(reshape2)
library(Matrix)
library(parallel)
library(future)
library(future.apply)
library(tidyverse)
library(nichenetr)
# library(Seurat)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(scales)
library(rlist)
library(stringr)
library(ggnet)
library(network)
library(ComplexHeatmap)
source("/public/home/gzcui_gdl/slurm.script/cellchat.visualization.R")
Ncore <- detectCores()
options(stringsAsFactors = FALSE)

packageVersion("CellChat")  ##‘1.5.0’

#### input nichenet database---
ligand_target_matrix <- readRDS("/public/home/gzcui_gdl/refgenome/nichenet/ligand_target_matrix.rds")
# str(ligand_target_matrix)
head(ligand_target_matrix[,1:4])

## input regulon information
regls <- readRDS("/public/home/gzcui_gdl/LCM/organogenesis/scenic/pyS6/regulon_stage.auc/E3.5_E8.75_s17.regulon.stage_geneTarget.ls.4678.rds")
head(str(regls))

regm <- "/public/home/gzcui_gdl/LCM/organogenesis/scenic/pyS6/regulon_stage.auc/regulon_stage.s17_MGI.phenotype_AUC_4678_3544s.auc.rds"
regm <- as.data.frame(readRDS(regm)) 
head(regm[,1:4])


# fls <-  list.files( path = paste0("."), pattern ="TPM.STGname_sampleFilt.rds",full.names = T) 
# fls

rdata <-  readRDS("/public/home/gzcui_gdl/LCM/organogenesis/scenic/matrix.all.stage/E3.5_E8.75_S17_3544sample_52812gene_TPM.STGname_sampleFilt.rds")
# str(rdata)
# rdata %>% dplyr::glimpse()
head(rdata[,1:4])


degsf <- readRDS("/public/home/gzcui_gdl/LCM/organogenesis/matrix.pre/deg.meg/LCM.E5.5_E8.75.domain.DEGs.filter.JUSTpVal001.316168.rds")

pho <- read.csv("/public/home/gzcui_gdl/LCM/organogenesis/trajectory/E3.5__E8.75_S17_3606s.domain_New.color.metaAll.scadj.csv",row.names = 1)
head(pho)
length(unique(pho$domain))
# setdiff(unique(pho$domain),unique(degsf$domain))
setdiff(unique(pho$orig.ident),NA)
# stages="E6.25"


##############################################################
ff <- setdiff(unique(pho$orig.ident),NA)[grep(stages,setdiff(unique(pho$orig.ident),NA))]

# ff <- "E5.5_181126_TPMlog2"
ffs <- gsub("_TPMlog2","",ff)
print(paste0("The stage running is ",ffs))

fn <- paste0("cellchat.Domain.log2TPM.",ffs)

# dir.create(ffs)
setwd(ffs)

pho1 <- filter(pho,orig.ident==ff)
str(pho1)
rda1 <- rdata[, colnames(rdata) %in% pho1$stname]
str(rda1)
stg <- unique(pho1$stage)



### run liana for enrich L-Rs---------
if (all(rownames(pho1)%in% colnames(rda1) ) & all(colnames(rda1) %in% rownames(pho1)) ) {
  print(paste0("The metadata is mattched with exp matrix in",ffs))
  
} else {
  print(paste0("The metadata is not mattched with exp matrix in",ffs))
  pho1 <- pho1[rownames(pho1) %in% colnames(rda1),]
  
}



data.input = log2(rda1+1) ## log(TPM) as input
fil <- apply(data.input, 1, function(x) length(x[x>1])>=5)
summary(fil)
# which(fil==FALSE)
data.input<-as.matrix(data.input[fil,])
head(data.input[,1:4])
data.input=na.omit(data.input)

#normalized
# data.input <-  normalizeData(data.input,scale.factor = 10000, do.log = TRUE) ##
str(data.input)
range(data.input)

meta <- data.frame(labels=pho1$domain,sample=pho1$stname,row.names = pho1$stname)
meta <- meta[colnames(data.input),]
head(meta)
# 


####**********************************************************************************
#### Step 1: CellChat L-R Pairs Loading
#### CellChat previously calculated enriched ligand-receptor (L-R) pairs between 
#### source and target domains. Here we load these results for downstream NicheNet analysis.
####**********************************************************************************

####------------------------------------------------
## NicheNet for downstream analysis
####------------------------------------------------
# Load CellChat output containing enriched L-R pairs between domains
fls <-  list.files( path = paste0("."), pattern ="df.net.all.csv",full.names = T) 
print(fls)
df.net <- read.csv(fls,row.names = 1)
head(df.net)

# Extract unique ligands from CellChat results and map to NicheNet database
ligands <- toupper(unique(df.net$ligand))  
length(ligands)
ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
length(ligands) 

# Filter L-R pairs to only include ligands present in NicheNet database
# and exclude autocrine interactions (source != target)
dfnich <- filter(df.net, toupper(ligand) %in% ligands )
dfnich <- filter(dfnich, source != target)
dim(dfnich)

# Create unique source-target pair identifiers for iteration
dfnich$pairs <- paste0(dfnich$source,"__",dfnich$target)
print( unique(dfnich$pairs))
# length(unique(dfnich$ligand))

#### Step 2: Iterate over each source-target domain pair for NicheNet analysis
for (p in unique(dfnich$pairs)) {
  # p="E6.25_exEpi__E6.25_Epi" 
  dir.create(p)
  
  # Extract L-R pairs for current source-target pair
  nich1 <- filter(dfnich, pairs == p)
  str(nich1)
  
  str(data.input)
  colnames(data.input)
  str(meta)
  
  #### Step 2.1: Define background_expressed_genes for NicheNet
  #### Genes expressed in target domains with log2(TPM+1) > 1 and expressed in at least 3 samples
  fil<- data.input[, rownames(meta)[meta$labels == as.character(unique(nich1$target)) ] ] %>%
    apply( 1, function(x) length(x[x>1])>=3) 
  background_genes <-  rownames(data.input)[fil]
  # length(background_genes)
  print(paste("The background genes in ",as.character(unique(nich1$target)), "are as follow:" )) 
  str(background_genes)
  
  #### Step 2.2: Define gene set of interest (geneset_oi) - Differentially Expressed Genes (DEGs) in target domains
  # geneDEG <- filter(degsf,domain ==  "E5.5_Epi" )
  geneDEG <- filter(degsf,domain ==  as.character(unique(nich1$target)) )
  geneDEG <- geneDEG[order(geneDEG$avg_log2FC,decreasing = T),]
  print(head(geneDEG)) 
  geneDEG$gene[1:30]
  # topdegnum=1000
  
  # Process DEGs only if available for the target domain
  if ( dim(geneDEG)[1]>0 ) {
    # Select top DEGs (up to topdegnum) and map to NicheNet database
    if (dim(geneDEG)[1] < topdegnum ) {
      
      geneset_oi <- geneDEG$gene %>%
        toupper() %>% 
        .[. %in% rownames(ligand_target_matrix)]
      
    } else {
      
      geneset_oi <- geneDEG$gene[1:topdegnum] %>%
        toupper() %>% 
        .[. %in% rownames(ligand_target_matrix)]
      
    }
    
    
    fn <- paste0("nichenet.cellchat.Domain.log2TPM_",ffs,"_pairs.",p,"_topDEG",length(geneset_oi))
    print(fn)
    # head(geneset_oi)
    # length(geneset_oi)
    
    #### Step 2.3: NicheNet ligand activity analysis
    #### predict_ligand_activities function performs ligand activity analysis to assess
    #### which ligands have the highest potential to regulate the gene set of interest (DEGs)
    ligand_activities <- predict_ligand_activities(
      geneset = geneset_oi,                    # DEGs in target domain (gene set of interest)
      background_expressed_genes = background_genes,  # Expressed genes in target domain
      ligand_target_matrix = ligand_target_matrix,    # NicheNet ligand-target prior model
      potential_ligands = unique(toupper(nich1$ligand))  # Ligands from CellChat L-R pairs
    )
    
    str(ligand_activities)
    saveRDS(ligand_activities,file = paste0(p,"/ligand_activities.nichenat_",fn,".rds"))
    # summary(ligand_activities)
    
    
    # ##-- select top liginds--------------
    # ligand_activities <- nichenet_activities
    # ligand_activities %>% arrange(-aupr_corrected)
    # 
    # if (dim(ligand_activities)[1]>90 ) {
    #   ntop <- round(dim(ligand_activities)[1]/3)
    # } else {
    #   ntop <- round(dim(ligand_activities)[1]/2)
    # }
    # print(paste0("The number of top ligand activities is :",ntop))
    # 
    # best_upstream_ligands_pn = ligand_activities %>% top_n(ntop, abs(aupr_corrected) ) %>% arrange(-aupr_corrected)
    # # best_upstream_ligands_pn = ligand_activities
    # head(best_upstream_ligands_pn)
    # best_upstream_ligands <- best_upstream_ligands_pn$test_ligand
    # str(best_upstream_ligands)
    
    # str(best_upstream_ligands_pn)
    # summary(best_upstream_ligands_pn)
    
    # show histogram of ligand activity scores
    # if (best_upstream_ligands_pn$aupr_corrected[ntop]<0) {
    #  pline <- min(best_upstream_ligands_pn$aupr_corrected[best_upstream_ligands_pn$aupr_corrected > 0] )
    #  nline <- max(best_upstream_ligands_pn$aupr_corrected[best_upstream_ligands_pn$aupr_corrected < 0] )
    # 
    #  p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    #    geom_histogram(color="black", fill="darkorange")  + 
    #    # geom_density(alpha=.1, fill="orange") +
    #    geom_vline(aes(xintercept=pline), color="red", linetype="dashed", size=1) + 
    #    geom_vline(aes(xintercept=nline), color="darkblue", linetype="dashed", size=1) + 
    #    labs(x="ligand activity (PCC)", y = "# ligands",title = paste0("Top_lig_activities:",ntop)) +
    #    theme_classic()   
    # } else {
    #   p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    #     geom_histogram(color="black", fill="darkorange")  + 
    #     # geom_density(alpha=.1, fill="orange") +
    #     geom_vline(aes(xintercept=min(ligand_activities %>% top_n(ntop, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
    #     labs(x="ligand activity (PCC)", y = "# ligands",title = paste0("Top_lig_activities:",ntop)) +
    #     theme_classic()
    # }
    # 
    # pdf(file = paste0(p,"/hist_lig_activity.",fn,".pdf"))
    # print(p_hist_lig_activity) 
    # dev.off()
    
    ## Use all ligands from ligand activity analysis
    best_upstream_ligands_pn = ligand_activities
    head(best_upstream_ligands_pn)
    best_upstream_ligands <- best_upstream_ligands_pn$test_ligand
    print(paste("use all ligands for nichenet analysis, the ligands are as follow:"))
    str(best_upstream_ligands)
    
    
    #### Step 2.4: Infer ligand-target regulatory links
    #### get_weighted_ligand_target_links function (n = 250) infers which genes in the gene set of interest
    #### have the highest regulatory potential for each top-ranked ligand
    active_ligand_target_links_df = best_upstream_ligands %>% 
      lapply(get_weighted_ligand_target_links, 
             geneset = geneset_oi,              # DEGs in target domain
             ligand_target_matrix = ligand_target_matrix,  # NicheNet ligand-target prior model
             n = 250) %>%                       # Top 250 target genes per ligand with highest regulatory potential
      bind_rows() %>% na.omit()
    print(paste("the active_ligand_target_links_df are: ",nrow(active_ligand_target_links_df))) 
    ## [1] 511 /1055
    head(active_ligand_target_links_df)
    # as.data.frame(active_ligand_target_links_df)
    
    # Prepare ligand-target matrix for visualization with cutoff = 0.25
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
    print(paste("the active_ligand_target_links are: ",nrow(active_ligand_target_links))) 
    ## [1] 91 /113
    head(active_ligand_target_links)
    str(active_ligand_target_links)
    saveRDS(active_ligand_target_links,file = paste0(p,"/active_ligand_target_links.nichenat_",fn,".rds"))
    
    if (nrow(active_ligand_target_links) >0) {
      
      order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
      order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(., rownames(active_ligand_target_links))
      
      
      vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
      str(vis_ligand_target)
      range(vis_ligand_target)
      saveRDS(vis_ligand_target,file = paste0(p,"/res1.heatmap_vis_ligand_target.nichenat_",p,".rds"))
      
      #####------plot in heatmap-------------
      # fn <- paste0(fn,"2")
      ytext <- paste0("Source_ligands:",as.character(unique(nich1$source)))
      xtext <- paste0("Signature in ",as.character(unique(nich1$target)))
      
      nich2 <- nich1[,c("ligand","pathway_name")]
      nich2$ligand <- toupper(nich2$ligand)
      nich2 <- unique(nich2)
      ylabs <- nich2$pathway_name
      names(ylabs) <- nich2$ligand
      str(ylabs)
      
      p_ligand_target_network = vis_ligand_target %>%
        make_heatmap_ggplot( ytext, xtext, 
                             color = "purple",legend_position = "right", 
                             x_axis_position = "top",legend_title = "Regulatory potential") +
        scale_y_discrete( labels =ylabs)+
        scale_fill_gradient2(low = "lightgreen",  high = "purple",mid="whitesmoke") + theme(axis.text.x = element_text(face = "italic"))
      
      # pdf(file = paste0(p,"/heatmap_ligand_target_network.",fn,".pdf"),width = round(dim(vis_ligand_target)[1]/2) ,height = round(dim(vis_ligand_target)[1]/4))
      # print(p_ligand_target_network)
      # dev.off()
      
      ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
      vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
      # str(vis_ligand_aupr)
      # head(ligand_aupr_matrix)
      
      p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot(ytext,"Ligand activity", 
                                                              color = "darkorange",
                                                              legend_position = "right", x_axis_position = "top",
                                                              legend_title = "AUPR\n(target gene prediction ability)")+
        scale_fill_gradient2(low = "lightgreen",  high = "darkorange",mid="whitesmoke")
      
      
      # pdf(file = paste0(p,"/heatmap_ligand_aupr.",fn,".pdf"),width =5 ,height = round(dim(vis_ligand_aupr)[1]/4))
      # print(p_ligand_aupr)
      # dev.off()
      
      # Combine the different heatmaps in one overview figure
      figures_without_legend = plot_grid(
        p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
        # p_ligand_tumor_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
        p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""), 
        align = "hv",
        nrow = 1,
        rel_widths = c(ncol(vis_ligand_aupr)+ 4.5,  ncol(vis_ligand_target)-2) 
      ) 
      
      legends = plot_grid(
        as_ggplot(get_legend(p_ligand_aupr)),
        as_ggplot(get_legend(p_ligand_target_network)),
        nrow = 2,
        align = "v")
      
      pdf(file = paste0(p,"/combine.heatmap_ligand_target_network.",fn,".pdf"),width = max(round(dim(vis_ligand_target)[1]),8),height = max(round(dim(vis_ligand_target)[1]/4),4) )
      p_plot_grid <- plot_grid(figures_without_legend, 
                               NULL,
                               legends, 
                               rel_widths = c(20,1,1), nrow = 1,ncol=3, align = "h")
      print(p_plot_grid)
      dev.off()
      #### end plot---------
      
      ####********************************************************************************
      #### Step 3: SCENIC Regulon Integration
      #### SCENIC reconstructs regulons (transcription factors and their target genes).
      #### Here we connect NicheNet-predicted target genes to their upstream regulons,
      #### and calculate the proportion of domain-specifically highly expressed target genes 
      #### among all target genes of each regulon.
      ####********************************************************************************
      
      ### ligand-receptor-target-regulon pathway construction------------------------------
      str(active_ligand_target_links) ## row- target; col-ligand
      str(best_upstream_ligands)
      str(nich1)
      length(unique(nich1$ligand) )
      length(unique(nich1$receptor) )
      length(unique(rownames(active_ligand_target_links)))
      targ <- unique(rownames(active_ligand_target_links))  # Target genes from NicheNet
      
      ## Select stage-specific regulons that contain NicheNet-predicted targets
      reg.sg <- regls[grepl(stg,names(regls))] ### E5.25_E5.5, E5.5_E5.75 regulons
      summary(grepl(targ,toupper(reg.sg) ))
      head(toupper(reg.sg[1]))
      
      # Extract regulons that contain at least one NicheNet target gene
      reg.taronly <- list()
      j=1
      for (i in seq_len(length(reg.sg)) ) {
        if (length(intersect( toupper( reg.sg[[i]]),targ))>0) {
          reg.taronly[[j]] <-  intersect( toupper( reg.sg[[i]]),targ)
          names(reg.taronly)[j] <- names(reg.sg[i])
          j=j+1
        } 
      }
      print(" The regulons enriched with NicheNet targets are as follow: ")
      str(reg.taronly)
      
      reg.tar <- list.subset(reg.sg, names(reg.sg) %in% names(reg.taronly))
      str(reg.tar)
      
      regname <- gsub("__.*$","",names(reg.tar))
      regname <- unique(regname)
      head(regname)
      length(regname)
      
      #### Calculate the proportion of domain-specific target genes among all target genes of each regulon
      #### This measures how enriched each regulon is with NicheNet-predicted targets
      tar.ratio <- lengths(reg.taronly)/lengths(reg.tar)
      head( sort(tar.ratio,decreasing = T) )
      summary(tar.ratio)
      # summary(lengths(reg.taronly))
      
      ## Calculate average AUC scores for regulons in the target domain
      ## AUC (Area Under Curve) represents regulon activity (transcription factor activity)
      head(pho1)
      tardomain <- unlist( strsplit(p,"__"))[2]
      tar.sam <-pho1  %>%  filter(domain==tardomain) %>% rownames()
      
      reg.auc <- regm[names(reg.tar),tar.sam]
      # head(reg.auc[,1:4])
      # dim(reg.auc)
      ave.auc <- apply(reg.auc, 1, mean)  # Average regulon activity across target domain samples
      head(sort(ave.auc,decreasing = T))
      # ave.auc["Ahr__E5.25_E5.5"]
      
      # Combine target ratio and average AUC for each regulon
      reg.tar.ratio.auc <- data.frame(row.names = names(tar.ratio),ratio.tar=tar.ratio)
      head(reg.tar.ratio.auc)
      reg.tar.ratio.auc <- cbind(reg.tar.ratio.auc,ave.auc)
      head(reg.tar.ratio.auc)
      
      
      saveRDS(reg.tar,file = paste0(p,"/nichenet.target.All.regulonList.rds") )
      saveRDS(reg.taronly,file = paste0(p,"/nichenet.target_ONLY.regulonList.rds") )
      
      #### Step 4: Merge ligand-receptor-target-regulon into complete signaling pathways
      ### Combine CellChat L-R pairs, NicheNet ligand-target links, and SCENIC regulon information
      lig.tar <- melt(active_ligand_target_links,value.name="prob.lig_tar")
      colnames(lig.tar)[c(1:2)] <- c("targetGene","ligand")
      lig.tar <- filter(lig.tar,prob.lig_tar>0)
      lig.tar$ligand <-  str_to_title(lig.tar$ligand)
      lig.tar$targetGene <-  str_to_title(lig.tar$targetGene)
      dim(lig.tar)
      # length(unique(lig.tar[,1]))
      head(lig.tar)
      # summary(lig.tar$prob.lig_tar)
      
      sigflow <- merge(nich1,lig.tar,by="ligand")
      summary(sigflow$prob.lig_tar)
      head(sigflow)
      dim(sigflow)
      saveRDS(sigflow,file = paste0(p,"/sigflow.all.rds") )
      
      reg.tar.pair <- data.frame(targetGene=NULL,regulon=NULL)
      for (i in seq_len(length(reg.taronly)) ) {
        # i=1
        pa <- data.frame(targetGene=reg.taronly[[i]],regulon=names(reg.taronly)[i])
        # dim(pa)
        # head(pa)
        reg.tar.pair <- rbind(reg.tar.pair,pa)
        
      }
      
      
      reg.tar.ratio.auc$regulon <- rownames(reg.tar.ratio.auc)
      reg.tar.pair <- merge(reg.tar.pair,reg.tar.ratio.auc,by="regulon")
      reg.tar.pair$targetGene <-  str_to_title( reg.tar.pair$targetGene)
      dim(reg.tar.pair)
      head(reg.tar.pair)
      write.csv(reg.tar.pair,file = paste0(p,"/target.regulons.",tardomain,"ratio.aveAuc.csv"))
      
      
      reg.tar.pair.f <- filter(reg.tar.pair,ratio.tar > summary(reg.tar.pair$ratio.tar)[5])
      dim(reg.tar.pair.f)
      
      sigflow <- merge(sigflow,reg.tar.pair.f,by="targetGene",all.y=T)
      head(sigflow)
      dim(sigflow)
      saveRDS(sigflow,file = paste0(p,"/sigflow.regulon.filter",p,".rds") )
      
      cat(paste(
        "the number of ligand :", length(unique(sigflow$ligand)), "\t",
        "the number of receptor :", length(unique(sigflow$receptor)),"\t",
        "the number of targets :",  length(unique(sigflow$targetGene)),"\t",
        "the number of regulon :", length(unique(sigflow$regulon))
      ),file = paste0(p,"/signaling.flow.information.txt"),append = T,sep = "\t")
      
      
      nlig <- length(unique(sigflow$ligand))
      nrep <- length(unique(sigflow$receptor))
      ntar <- length(unique(sigflow$targetGene))
      nreg <- length(unique(sigflow$regulon))
      ny.max <-  max(nlig,nrep,ntar,nreg)
      
      ####################################################################################
      #### Step 5: Signaling Network Construction and Visualization
      #### Construct a 4-layer signaling network: Ligand -> Receptor -> Target -> Regulon
      ####################################################################################
      
      head(sigflow)
      
      # Create network edges for each layer:
      # Layer 1: Ligand (L:) -> Receptor (R:) with CellChat probability weight
      edgs1 <- unique(data.frame(from=paste0("L:", sigflow$ligand),to= paste0("R:",str_to_title(sigflow$receptor)),weight=sigflow$prob,lty=1)) 
      
      # Layer 2: Receptor (R:) -> Target gene (T:) with NicheNet ligand-target probability
      edgs2 <- unique(data.frame(from=paste0("R:",str_to_title(sigflow$receptor)),to= paste0("T:",sigflow$targetGene) ,weight=sigflow$prob.lig_tar,lty=3)) 
      
      # Layer 3: Regulon (reg:) -> Target gene (T:) with average weight from previous layers
      # This represents transcription factor regulation of target genes
      edgs3 <- unique(data.frame(from=paste0("reg:",sigflow$regulon) ,to=  paste0("T:",sigflow$targetGene),weight=mean( mean(edgs1$weight), mean(edgs2$weight)),lty=1)) 
      
      # Combine all edges
      edgs <- rbind(edgs1,edgs2,edgs3)
      
      # Remove duplicate edges
      dul <- duplicated(paste0(edgs$from,edgs$to))
      summary(dul)
      # edgs[edgs$from =="R:Bmpr1a_bmpr2"& edgs$to=="T:Zfhx3" ,]
      edgs <- edgs[!dul,]
      dim(edgs)
      
      
      # Position nodes in 4 layers (columns) for clear visualization:
      # Layer 1 (X=1): Ligands (colored coral)
      node1 <-  arrange(edgs1, desc(weight)) 
      ny <- unique(node1$from)
      nstar <- (ny.max-length(ny))/2
      nd1 <- data.frame(node=ny,X=1,Y= seq(-nstar,(-length(ny)-nstar+1),by=-1),colr="coral",size=1)
      
      # Layer 2 (X=2): Receptors (colored burlywood1)
      ny <- unique(node1$to)
      nstar <- (ny.max-length(ny))/2
      nd2 <- data.frame(node=ny,X=2,Y=seq(-nstar,-length(ny)-nstar+1,by=-1),colr="burlywood1",size=1)
      
      # Layer 3 (X=3.5): Target genes (colored darkolivegreen3)
      node2 <-  arrange(edgs2, desc(weight)) 
      ny <-unique(node2$to)
      nstar <- (ny.max-length(ny))/2
      nd3 <- data.frame(node=ny,X=3.5,Y=seq(-nstar,-length(ny)-nstar+1,by=-1),colr="darkolivegreen3",size=1)
      
      # Layer 4 (X=5): Regulons (colored by target ratio, sized by AUC activity)
      node3 <- unique(data.frame(reg=sigflow$regulon, auc=sigflow$ave.auc,ratio=sigflow$ratio.tar)) 
      node3 <-  arrange(node3, desc(auc)) 
      node3$sal.auc <- round(rescale(node3$ratio,to=c(1,200))) 
      clrs <- colorRampPalette(c("plum1", "purple"))(200)
      node3$colr <- clrs[node3$sal.auc]
      ny <- paste0("reg:",unique(node3$reg)) 
      nstar <- (ny.max-length(ny))/2
      nd4 <- data.frame(node=ny,X=5,Y=seq(-nstar,-length(ny)-nstar+1,by=-1),colr=node3$colr,size=node3$auc)
      
      node.xy <- rbind(nd1,nd2,nd3,nd4)
      rownames(node.xy) <- paste0(node.xy$node)  
      head(node.xy)
      dim(node.xy)
      
      all(node.xy$node %in% unique(c(edgs$from,edgs$to)))
      
      cord.xy <- node.xy[,c("node","X","Y")]
      rownames(cord.xy) <- cord.xy$node
      cord.xy <- as.matrix(cord.xy[,-1])
      head(cord.xy)
      
      
      net <- as.network(edgs, directed = T, vertices = node.xy,multiple = F )
      # net %v% "exp.color" <- brewer.pal(9 , 'Purples')[net %v% "Exp.ave" ]
      set.edge.attribute(net, "edg.color", "gray30" )
      # set.edge.attribute(net, "edg.color", ifelse(net %e% "pos.neg" %in% c("UP","Down"), ifelse( net %e% "pos.neg" =="Down", "green4","red4"), "lightblue" ))
      
      max(nlig,nrep,ntar,nreg)
      saveRDS(net, file = paste0(p,"/net.project.",p,".rds") )
      
      pp <- ggnet2(net,mode = cord.xy,
                   color="colr", shape=19, #color="Exp.ave",palette="Purples"
                   size="size", max_size=15,
                   size.legend="AUC",color.legend="type.enrich",legend.position="none",
                   label=T,label.size = 5,
                   edge.size = "weight",edge.color = "edg.color",edge.alpha=0.6, edge.lty = "lty",
                   arrow.size = 2, arrow.gap = 0.2)
      
      pdf(file = paste0(p,"/signal.flow.network", fn,"net.pdf"),width = 10,height=max(length(unique(sigflow$ligand)),length(unique(sigflow$receptor)),length(unique(sigflow$targetGene)),length(unique(sigflow$regulon)))/2,5 )
      print(pp)
      dev.off() 
      
    }
    
    
    
  } else {
    print( paste("There are no DEGs in domain: ",as.character(unique(nich1$target))," Skip the nichnet analysis" )  )
  }
  
}




sessionInfo()

q()















