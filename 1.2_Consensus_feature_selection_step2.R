##Geo-seq data, script2--Find the top PC30 feature names from expression matrix
##Guizhong Cui
## 30-9-22
# rm(list=ls())

args <- commandArgs(trailingOnly = T)

ff <- args[1]
dir <- args[2] #workdir qsub E.. 
nfac <-as.integer(args[3])  #the number of NMF factors
npc <- as.integer(args[4]) #the number of PC dims
kmin <-as.integer(args[5])  #the min number of KNN 
kmax <-as.integer(args[6]) 

print(ff)

library(future)
library(future.apply)
library(parallel)
plan("multicore", workers = detectCores())
print(paste0("the number of core used ",detectCores()))

# dim(Matrix)
fls <-  list.files(path =dir, pattern ="samples_Newname_geneFilt.sampleFilt.rds",
                   recursive = TRUE,full.names = F) 
Matrix <- readRDS(paste0(dir,"/",fls))
dim(Matrix)
range(Matrix)

fls <-  list.files(path =dir, pattern ="samples_Newname_geneFilt.sampleFilt_regress.ngene.pmt.noscale.rds",
                   recursive = TRUE,full.names = F) 
da.reg <- readRDS(paste0(dir,"/",fls))
dim(da.reg)
range(da.reg)
head(da.reg[,1:5])

####====get duplicated features================================================================
fn <- paste0(ff,"_5FS_")

library(VennDiagram)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
source("../help_code/Plot.funtions.lcm.R")

# dir.create("3.skfeature")
#linux run:
# mv fs_* 3.skfeature/
# mv geneName_fs_* 3.skfeature/


files<-dir("3.skfeature",pattern = "geneName_fs_")
fs.s <- list()
for ( i in seq_len(length(files))){
  # i=1
  tempcsv<-read.csv(paste0("3.skfeature/",files[i]),header = T,row.names = 1)
  head(tempcsv)
  fs.s[i] <- list(unique(tempcsv[,1])) 
  fs.mth <- gsub(paste0("_",ff,".*$"),"",files[i])
  # fs.mth <- gsub(".2500_.*$","",files[i])
  # fs.mth <- gsub(".3000_.*$","",files[i])
  names(fs.s)[i] <-gsub("geneName_fs_","",fs.mth)
  # str(fs.s)
  
}
str(fs.s)

p.venn(degs = fs.s, fn= paste0("3.skfeature/Venn_multi.featureSlection_"))

fs.ls <- as.character(unlist(fs.s[c(1:5)]))  
head(fs.ls) 
# str(fs.ls)
fs.du <- unique(fs.ls[duplicated(fs.ls)] ) 
length(fs.du) ##3104
sample(fs.du, 20)
# fs.u <- unique(fs.ls)
# length(fs.u) #5720
dim(Matrix)
fsdu <- fs.du[(fs.du %in% rownames(da.reg))]
length(fsdu) #3059

setdiff(fs.du,fsdu)

# fs.matrix <- Matrix[fsdu,]
fs.matrix <- da.reg[fsdu,]
# dim(fs.matrix)

if (sum(which(is.na(fs.matrix)))==0 & all(rownames(fs.matrix %in% fsdu))) {
  saveRDS(fs.matrix,file = paste0("PC30selected_5FSduplicated_",dim(fs.matrix)[1],"genes_",fn,".rds"))
  saveRDS(fsdu,file = paste0("PC30selected_5FSduplicated_",dim(fs.matrix)[1],"geneslist.",fn,".rds"))
} else { print("NA in the fs.matrix")}



#######--NMF feature extraction---------------------------------------------------------
fn <- paste0(ff,"_FS_nmfSelect_")
# library(CelliD)
library(tidyverse) # general purpose library for data handling
library(ggpubr) #library for plotting
library(leiden)
library(clustree)
library(Seurat)
packageVersion("Seurat") #'4.1.1'
packageVersion("clustree") #0.5.0
packageVersion("dplyr") #1.1.2'

##meta information
fmeta <-  list.files(path =dir, pattern ="filter.meta.csv",
                     recursive = TRUE,full.names = F) 

Meta <- read.csv(paste0(fmeta),row.names = 1)
head(Meta)
dim(Meta)
pho <- Meta[Meta$filt.all==TRUE,c("newname","seqname","section","ngene","percent.mt","sector23","lineage","colr.s","colr.l")]
str(pho)

## color setting---
cols1 <-unique( pho$colr.l)
names(cols1) <- unique(pho$lineage)

cols2 <- unique(pho$colr.s)
names(cols2) <- unique(pho$sector23)


#constructed object, all expressed gene-------------------
pbmc <- CreateSeuratObject(counts = Matrix,min.cells = 2, meta.data =pho)
pbmc
# pbmc <- CreateSeuratObject(counts = da.reg,min.cells = 2, meta.data =pho)
# pbmc

dir.create("4.cluster")
fsgn <- fsdu
fsgn <- fsgn[(fsgn %in% rownames(pbmc))]
str(fsgn)
# Library-size normalization, log-transformation, and centering and scaling of gene expression values
# pbmc<- NormalizeData(pbmc)  #log=TRUE was slected
pbmc <- ScaleData(pbmc, vars.to.regress = c("ngene", "percent.mt"), features = rownames(pbmc))
# pbmc <- ScaleData(pbmc, features = rownames(pbmc))
# pbmc <- ScaleData(pbmc, features = fsgn)
str(pbmc)


###----NMF feature extraction--------------------------------------------------------
# nfac <- 11
alldim <-nfac  # the number of NMF factors
pbmc<- STutility::RunNMF(pbmc,nfactors =alldim,features = fsgn)

# str(pbmc)
####----factors selection ---NMF----------------------
# str(pbmc)
nmftop <- pbmc@reductions$NMF@cell.embeddings
str(nmftop)
# head(nmftop[,1:5])
# summary(nmftop)

##--density of NMF factor
pdf(file = paste0("4.cluster/","nmf_sample_density.Factor",alldim,"_",fn,".pdf"))
for (i in 1:alldim) {
  # i=1
  fac <- nmftop[,i]
  d <- density(fac)
  plot(d, main=paste0("nmf_density.Factor",i))
}
dev.off()

#--kemean=2 of NMF factor--
fnmf <- list()
for (i in 1:alldim) {
  # i=1
  fac <- nmftop[,i]
  km <- kmeans(fac, 2, nstart = 1)
  kmc <- km$cluster
  
  if (mean(fac[kmc==1])>mean(fac[kmc==2])) {
    snmf <- list(names(fac)[kmc==1]) 
    fnmf[i] <- snmf
  } else {
    fnmf[i] <- list( names(fac)[kmc==2])
  }
  names(fnmf)[i] <- i
  
}
str(fnmf)
fl.nmf <- lengths(fnmf)
nmf.fc<- fl.nmf[fl.nmf>=5]
nmf.dim <- as.numeric(names(nmf.fc)) 
length(nmf.dim)

nmfinfo <- as.data.frame(summary(fnmf)) 
cat(paste0("Total.NMF.Factors:",length(nmf.dim),"; \t"),file = paste0("4.cluster/nmf.selected.factor.info.",fn,".parameter.txt"),append = T,sep = "\t")
cat(paste0("sample number of each nmf factors: \t",nmfinfo, "\t") ,file = paste0("4.cluster/nmf.selected.factor.info.",fn,".parameter.txt"),append = T,sep = "\t")


####----gene selection ---NMF----------------------
# str(pbmc)
nmftop <- pbmc@reductions$NMF@feature.loadings
str(nmftop)
head(nmftop[,1:5])
# summary(nmftop)

##--density of NMF factor
pdf(file = paste0("4.cluster/","nmf_gene_density.Factor",alldim,"_",fn,".pdf"))
for (i in nmf.dim) {
  # i=1
  fac <- nmftop[,i]
  d <- density(fac)
  plot(d, main=paste0("nmf_density.Factor",i))
}
dev.off()

#--kemean=2 of NMF factor--
# nmfls <- list()
# for (i in nmf.dim) {
#   # i=1
#   fac <- nmftop[,i]
#   km <- kmeans(fac, 2, nstart = 1)
#   kmc <- km$cluster
#   
#   if (mean(fac[kmc==1])>mean(fac[kmc==2])) {
#     snmf <- list(names(fac)[kmc==1]) 
#     
#     nmfls[i] <- snmf
#   } else {
#     nmfls[i] <- list( names(fac)[kmc==2])
#    }
#   names(nmfls)[i] <- paste0("NMF.factor",i)
#   
# }
# str(nmfls)
# 
# fs.nmf <- unique(unlist(nmfls))
# length(fs.nmf)


#--kemean=3 of NMF factor--
nmfls <- list()
for (i in nmf.dim) {
  # i=1
  fac <- nmftop[,i]
  km <- kmeans(fac, 3, nstart = 1)
  kmc <- km$cluster
  km.ave <- c(mean(fac[kmc==1]),mean(fac[kmc==2]),mean(fac[kmc==3]))
  kmthre <- (sort(km.ave)[2]+sort(km.ave)[3])/2
  
  if (mean(fac[kmc==1])>kmthre) {
    snmf <- list(names(fac)[kmc==1]) 
    
    nmfls[i] <- snmf
  } else if (mean(fac[kmc==2])>kmthre) {
    nmfls[i] <- list( names(fac)[kmc==2])
  } else {
    nmfls[i] <- list( names(fac)[kmc==3])
  }
  
  
  names(nmfls)[i] <- paste0("NMF.factor",i)
  
}
str(nmfls)

fs.nmf <- unique(unlist(nmfls))
nmfinfo <- as.data.frame(summary(nmfls)) 

cat(paste0("Total.NMF.Genes:",length(fs.nmf),"; \t"),file = paste0("4.cluster/nmf.selected.factor.info.",fn,".parameter.txt"),append = T,sep = "\t")
cat(paste0("gene number of each nmf factors: \t",nmfinfo, "\t") ,file = paste0("4.cluster/nmf.selected.factor.info.",fn,".parameter.txt"),append = T,sep = "\t")

saveRDS(fs.nmf,file = paste0("Featrue.extraction_NMF_basedPC30.5FSdup_NMFfactor",alldim,"_feature",length(fs.nmf),".rds"))
########-------------NMF feature extraction end--------------------------------------


###--pheatmap for fs.nmf-----------------------------------------------------------
library(pheatmap)
library(RColorBrewer)
fn1<-paste0("4.cluster/",ff,"_NMFfeature",length(fs.nmf),"_regress.heatmap")

#exp matrix for heatmap
head(da.reg[,1:4])
deg.m <- da.reg[fs.nmf,]
p1.d <- deg.m
dim(p1.d)
all(rownames(p1.d) %in% fs.nmf)
## meta info
# dim(meta)
anno <-pho[,c("sector23","lineage","section")]
# head(anno)
# summary(anno)
anno$section<-as.factor(anno$section)
# levels(anno$section)

# str(cls)
sec <- colorRampPalette(c('grey','black'))(nlevels(anno$section))
names(sec) <- levels(anno$section)
##anno color
ann_colors = list(sector23 = cols2,lineage=cols1,
                  # seurat_clusters=cls,
                  section=sec
                  
)

ng <- dim(p1.d)[1]
if (ng>2000) {
  chv <- 0.1
} else {
  if (ng >1000) {
    chv <- 0.2
  } else {
    if (ng >500) {
      chv <- 0.4
    } else {
      chv <- 1
    }
  }
}

ns <- dim(p1.d)[2]
if (ns>1000) {
  rhv <- 0.5
} else {
  if (ng >500) {
    rhv <- 1
  } else {
    if (ng >300) {
      rhv <- 1.5
    } else {
      rhv <- 2
    }
  }
}
print(c(ng,ns,fn1,chv,rhv))
p.heat(pda = p1.d,annotation = anno,color.anno = ann_colors,filename = fn1,ch = chv,rh=rhv)


###----PCA factors selection___parallel analysis--------------------
# library(psych)
# dim(pbmc@assays$RNA@scale.data[fs.nmf,])
# 
# pdf(file = paste0(fn1,"psych.PCparallel.pdf"))
# pa <- fa.parallel(pbmc@assays$RNA@scale.data[fs.nmf,], fa = "pc", n.iter = 100,
#                   show.legend = F, main = "Scree plot with parallel analysis")
# print(pa)
# dev.off()
# str(pa)
# saveRDS(pa, file = paste0(fn1,"psych.PCparallel.obj.rds"))
# 
# pdf(file = paste0(fn1,"psych.PCparallel.long.pdf"),width = 30,height = 5)
# plot(pa)
# dev.off()
# 
# pcs <- seq_len(pa$ncomp )
# cat(paste("parrallel analysis suggest top PC:",pcs,sep = " "),file = paste0("4.cluster/",fn,".parrallel.PC.selected.parameter.txt"),append = T)
############-----------Leiden pre-clustering-------------------------------------
fn <- paste0(ff,"_5FSnmf",length(fs.nmf),"_")
print(fn)

###dimensionality reduction ----------PCA-------
# pbmc<- RunICA(pbmc, features =  VariableFeatures(object = pbmc),nics = nics) ##ics 
# npc <- 30
pbmc<- RunPCA(pbmc, features = fs.nmf,npcs = npc)
str(pbmc)

pdf(file=paste0("4.cluster/",fn,"NMF",alldim,".factor.enrichmentScore.PCA.pdf"),width = 16,height = 12)
p <- FeaturePlot(pbmc, reduction = "pca", features =  paste0("factor_", 1:alldim), ncol = 4)
plot(p)
dev.off()


##significate dims--
pbmc <- JackStraw(object = pbmc, dims = npc,num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:npc)
# pbmc <- ScoreJackStraw(object = pbmc, dims = 21:40)
# str(pbmc)

pdf(paste0("4.cluster/",fn,".JackStrawPlot_",npc,"PCs.pdf"),width = 10,height = 6)
p <- JackStrawPlot(object = pbmc, dims = 1:npc)
plot(p)
dev.off()

pdf(paste0("4.cluster/",fn,".ElbowPlot_",npc,"PCs.pdf"))
p <- ElbowPlot(object = pbmc,ndims = npc)
plot(p)
dev.off()

##PC significance
pc.sig <- pbmc@reductions$pca@jackstraw@empirical.p.values #hvg*PC
# head(pc.sig[,1:4])
# head(Loadings(pbmc, reduction = "pca")[, 1:5])
pc.p <-as.data.frame(pbmc@reductions$pca@jackstraw@overall.p.values)  # p-value of PC
# head(pc.p)
# str(pc.p)
pc.s <- filter(pc.p,Score<0.01)

pcs <- pc.s$PC #PC selected by significance

npc <- length(pcs)

# pc contribution
pc.sd <-  Stdev(pbmc, reduction = "pca")
pc.perc <-sum(pc.sd[pcs])/sum(pc.sd)

cat(paste(npc,round(pc.perc,3),";PCs:",sep = "_"),file = paste0("4.cluster/",fn,".PC.selected.parameter.txt"),append = T)
cat(pcs,file = paste0("4.cluster/",fn,".PC.selected.parameter.txt"),append = T,sep = "\t")

###-----test K------
for (k in c(seq(kmin,kmax,1))) {
  pbmc <- FindNeighbors(object = pbmc, reduction = "pca",dims =pcs,k.param = k)
  # str(pbmc)
  
  ###Clustering---leiden------
  ##cluster numbers test with leiden
  pbmc <- FindClusters(object = pbmc, resolution = c(seq(0,2,.2)), graph.name="RNA_snn", algorithm = 4,random.seed = 0)
  head(pbmc@meta.data)
  # meta <-pbmc@meta.data
  # write.csv(meta,file = paste(fn,ng,"meta.csv"))
  # clut.res <- pbmc@meta.data[,c(paste0("RNA_snn_res.",seq(1,2,.2)))]
  # head(clut.res)
  # colnames(clut.res) <- gsub("RNA_snn_res.","K",colnames(clut.res))
  # head(clut.res)
  # str(clut.res)
  #cluster number in different resolutions
  pdf(paste0("4.cluster/",fn,".K",k,".clustertree.pdf"),width = 20,height = 30)
  p<-clustree(pbmc@meta.data, prefix = "RNA_snn_res.") +
    labs(title = paste(fn,".K",k))
  plot(p)
  dev.off()
  
  # pdf(paste0("4.cluster/",fn,".K",k,".clustertree.pdf"),width = 20,height = 30)
  # p<-clustree(clut.res, prefix ="K") +
  #   labs(title = paste(fn,".K",k))
  # plot(p)
  # dev.off()
  
  
  #cluster information in different resolutions
  
  pdf(paste0("4.cluster/",fn,".K",k,".clustertree.secotr.pdf"),width = 20,height = 30)
  p<-clustree(pbmc@meta.data, prefix = "RNA_snn_res.",node_label = "sector23",node_label_aggr = "label_position") +
    labs(title = paste(fn,".K",k))
  plot(p)
  dev.off()
}


saveRDS(pbmc,file = paste0("seurat.project.",fn,".rds"))
sessionInfo()
q()

####-----feature selection & extraction end----------------------------------------------------
