###############################################################################
# Trajectory Construction Using Differentially Expressed Regulons
# 
# This script constructs developmental trajectories across consecutive embryonic
# stages using SCENIC regulon activity scores. It integrates differential 
# regulon analysis, dimensionality reduction, KNN-based lineage inference,
# and Sankey visualization to model the progression of cell/domain states.
#
# Workflow Overview:
# 1. Load regulon AUC scores and metadata for each stage pair
# 2. Perform differential regulon analysis between domains using Seurat
# 3. Apply dimensionality reduction (PCA, UMAP, tSNE) using DE regulons
# 4. Use KNN-based lineage inference to predict ancestor-descendant relationships
# 5. Calculate median edge probabilities across bootstrap iterations
# 6. Generate interactive Sankey plots to visualize trajectories
#
# Input:
#   - Regulon AUC score matrices (aucScore.csv) for each stage pair
#   - Metadata file (domain.color.metaAll.csv) with domain annotations
#   - Command line: stage index (ks), working directory (dir)
#
# Output:
#   - Differential regulon tables and heatmaps
#   - PCA/UMAP/tSNE visualizations
#   - KNN lineage inference results (probability matrices)
#   - Edge probability tables for trajectory construction
#   - Interactive Sankey plots (HTML)
#

#
# Author: cgz
# Version: S17 (Differential Regulons)
###############################################################################

##########################
#### Step 1: Data Preparation and Differential Regulon Analysis ####
##########################

# Clean workspace
rm(list=ls())

# Load required libraries
library(Seurat)           # Single-cell analysis toolkit
library(future)           # Parallel processing
library(future.apply)     # Parallel apply functions
library(dplyr)            # Data manipulation
library(ggplot2)          # Visualization
library(Matrix)           # Sparse matrix operations
library(htmlwidgets)      # HTML widget saving
library(plotly)           # Interactive visualization (Sankey plots)
library(monocle3)         # Trajectory analysis (not actively used here)
library(ggpubr)           # Plot arrangement (ggarrange)
library(FNN)              # Fast nearest neighbor search (get.knnx)
library(reshape2)         # Data reshaping (melt)
library(scales)           # Scale functions
library(gplots)           # Plotting utilities
library(viridis)          # Color scales

# Configure parallel processing
plan("multicore", workers = detectCores())
print(paste0("Number of cores used: ", detectCores()))

# Source custom plotting functions
source("Plot.funtions.R")

time_point = paste0("E", c( seq(5.5, 8.75, 0.25)))
# time_point <- setdiff(time_point,"E5.75") 
# time_point <- setdiff(time_point,"E6.5") 
time_point[3] <- "E6.0"
time_point[7] <- "E7.0"
time_point[11] <- "E8.0"
time_point
time_point <- c("E3.5","E4.5","E5.25",time_point)
time_point
length(time_point)


args = commandArgs(trailingOnly=TRUE)

## ks is the stages from E7.5
ks = as.numeric(args[1])
dir <- args[2] 

## get merge meta files!!!!!!!!!!------------------
# mfls <- list.files(path =paste0(dir), pattern =paste0("domain.color.metaAll.csv"),
#                    recursive = TRUE,full.names = F)
mfls <- list.files(path ="../", pattern =paste0("domain.color.metaAll.csv"),
                   recursive = F,full.names = T)
# mfls
pho <- read.csv(mfls,row.names = 1)
head(pho)

# Filter out samples with NA cluster assignments
pho <- dplyr::filter(pho, cluster != "NA")
print(paste0("Metadata dimensions after filtering: ", dim(pho)[1], " samples x ", dim(pho)[2], " features"))


# Main loop: Process each consecutive stage pair
# Iterate through all stage transitions (E3.5->E4.5, E4.5->E5.25, ..., E8.5->E8.75)
for (kk in seq_len(length(time_point) - 1)) {
  
  print(paste0("Processing stage transition k = ", kk))
  
  # Define current stage pair (earlier stage = time_1, later stage = time_2)
  time_1 <- time_point[kk]
  time_2 <- time_point[kk + 1]
  print(paste0("Stage transition: ", time_1, " -> ", time_2))

  # Create output directory for this stage pair
  dir.create(paste0(time_1, "_", time_2), showWarnings = FALSE)

  # Base filename for outputs
  fn <- paste0(time_1, "_", time_2, "_regulon")

  #####---Load Regulon AUC Scores----------
  # Find AUC score file for this stage pair
  fls <- list.files(path = dir, pattern = paste0("aucScore.csv"),
                    recursive = TRUE, full.names = FALSE)
  fls <- fls[grep(paste0(time_1, "_", time_2), fls)]
  
  # Read regulon AUC scores (rows = samples, columns = regulons)
  aucell.out <- read.csv(fls, header = TRUE, check.names = FALSE, 
                         row.names = 1, sep = ",")
  
  # Clean regulon names by removing "(+)" suffix
  colnames(aucell.out) <- gsub("\\(\\+)$", "", colnames(aucell.out))

  # Transpose matrix: rows = regulons, columns = samples
  reg <- t(aucell.out)
  print("Regulon matrix preview:")
  print(head(reg[, 1:3]))
  print(paste0("Regulon matrix dimensions: ", dim(reg)[1], " regulons x ", dim(reg)[2], " samples"))

  #####---Filter Metadata for Stage Pair----------
  # Extract metadata for the two stages
  meta12 <- filter(pho, stage == time_1 | stage == time_2)
  print(paste0("Metadata dimensions for stages ", time_1, " and ", time_2, ": ", dim(meta12)[1], " samples"))

  # Verify sample matching between regulon matrix and metadata
  print(paste0("Samples in regulon but not in metadata: ", 
               length(setdiff(colnames(reg), meta12$stname))))
  print(paste0("Samples in metadata but not in regulon: ", 
               length(setdiff(meta12$stname, colnames(reg)))))

  # Align metadata with regulon samples
  if (length(setdiff(colnames(reg), meta12$stname)) == 0) {
    print("All regulon samples have metadata. Subsetting metadata...")
    meta.12r <- meta12[meta12$stname %in% colnames(reg), ]
    rownames(meta.12r) <- meta.12r$stname
  } else {
    # Some regulon samples missing metadata - log and filter
    missing_samples <- setdiff(colnames(reg), meta12$stname)
    print(paste0("Warning: ", length(missing_samples), " samples have no metadata"))
    cat(paste0(missing_samples, collapse = ","), 
        file = paste0(fn, "_sample.noMeta.txt"))
    
    # Filter both metadata and regulon matrix to matching samples
    meta.12r <- meta12[meta12$stname %in% colnames(reg), ]
    rownames(meta.12r) <- meta.12r$stname
    reg <- reg[, colnames(reg) %in% meta.12r$stname]
    print(paste0("Regulon matrix after filtering: ", dim(reg)[1], " x ", dim(reg)[2]))
  }

  #####---Special Handling for E7.0 Stage----------
  # E7.0 has additional batch/embryo information that needs to be incorporated
  if (time_2 == "E7.0" | time_1 == "E7.0") {
    print("Incorporating E7.0 embryo batch information")
    
    # Load E7.0 batch metadata
    batif <- read.csv("../E7.0_ref.meg.f.batch.csv", header = TRUE)
    colnames(batif)[2] <- "embryo"
    
    # Create sample name matching pattern
    batif$stname <- paste0("E7.0_", batif$newname)
    batif <- batif[, c("stname", "embryo")]
    
    # Merge batch info with metadata
    meta.12r <- merge(meta.12r, batif, by = "stname", all.x = TRUE)
    meta.12r[is.na(meta.12r)] <- "other"  # Fill missing embryo info
    rownames(meta.12r) <- meta.12r$stname
  } 

  print(paste0("Final metadata dimensions: ", dim(meta.12r)[1], " samples"))
  
  # Proceed only if metadata and regulon matrix have matching sample counts
  if (dim(meta.12r)[1] == dim(reg)[2]) {
    
    #####---Create Seurat Object----------
    # Construct Seurat object with regulon AUC scores as "counts"
    pbmc <- CreateSeuratObject(counts = reg, meta.data = meta.12r)
    
    #####---Data Scaling with Covariate Regression----------
    # Scale data and regress out unwanted variation
    # For E7.0 stages: regress out embryo batch effects
    # For other stages: regress out stage effects
    if (time_2 == "E7.0" | time_1 == "E7.0") {
      pbmc <- ScaleData(pbmc, features = rownames(pbmc), vars.to.regress = c("embryo"))
      fn <- paste0(time_1, "_", time_2, "_regulon.regressoutEmbryo.")
    } else {
      pbmc <- ScaleData(pbmc, features = rownames(pbmc), vars.to.regress = c("stage"))
      fn <- paste0(time_1, "_", time_2, "_regulon.regressoutStage.")
    }
    
    #####---Differential Regulon Analysis----------
    # Identify regulons that are differentially active between domains
    
    
    # Set identity class to 'domain' for comparison
    Idents(object = pbmc) <- 'domain'
        
    # Run FindAllMarkers to identify differentially active regulons
    # Test: Wilcoxon rank-sum test
    # Only consider regulons present in at least 50% of cells in at least one group
    test <- "wilcox"
    pbmc.markers <- FindAllMarkers(pbmc, test.use = test, slot = "data", 
                                   random.seed = 1, min.cells.group = 2,
                                   only.pos = FALSE, min.pct = 0.5, logfc.threshold = 0.01)
    
    print(paste0("Found ", nrow(pbmc.markers), " differential regulons"))
    
    #####---Save Differential Regulon Results----------
    # Save all differential regulons
    write.csv(pbmc.markers, 
              file = paste0(time_1, "_", time_2, "/", fn, "_", 
                            length(unique(pbmc.markers$gene)), "DER.",
                            round(length(unique(pbmc.markers$gene))/dim(pbmc)[1], 3), 
                            "rate.", test, ".csv"))
    
    # Extract top 10 regulons per domain by log2FC
    top10 <- pbmc.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    str(top10)
    write.csv(top10,file=paste0(time_1, "_", time_2,"/",fn,"_",dim(top10)[1],"Top10.DEGs.",test,".csv"))
    
    #####---Heatmap of Top Differential Regulons----------
    pdf(file = paste0(time_1, "_", time_2, "/", fn, ".DER.heatmap.top10.pdf"),
        height = length(levels(pbmc)) * 2, width = length(levels(pbmc)))
    p <- DoHeatmap(pbmc, features = top10$gene) +
      scale_fill_gradientn(colors = c("darkgreen", "white", "red"))
    plot(p)
    dev.off()
   
    #####---Filter Significant Differential Regulons----------
    # Filter regulons with adjusted p-value < 0.05
    deg.s <- filter(pbmc.markers, p_val_adj < 0.05)
    print(paste0("Significant differential regulons (padj < 0.05): ", nrow(deg.s)))
    
    # Extract scaled data (after regression) for significant regulons
    reg.rg <- GetAssayData(pbmc, slot = "scale.data")
    deg.w <- reg.rg[unique(deg.s$gene), ]
    
    # Log selection parameters
    cat(paste("Number of significant regulons (padj < 0.05):", 
              length(unique(deg.s$gene)), 
              "(", round(length(unique(deg.s$gene))/dim(pbmc)[1] * 100, 1), "% of all regulons)",
              sep = " "), 
        file = paste0(time_1, "_", time_2, "/", fn, ".regulon.selected.parameter.txt"),
        append = TRUE)
    
        #####---Dimensionality Reduction Parameters----------
    # Set stage-specific parameters for dimensionality reduction
    if (kk == 1) {
      # E3.5 -> E4.5 transition
      npc <- 15           # Number of principal components
      n.neighbors <- 3    # UMAP neighbors
      perplexity <- 10    # tSNE perplexity
    } else if (kk == 2) {
      # E4.5 -> E5.25 transition
      npc <- 15
      n.neighbors <- 10
      perplexity <- 10
    } else {
      # E5.25 -> later stages
      npc <- 30
      n.neighbors <- 20
      perplexity <- 15
    }
    
    #####---Select Variable Regulons----------
    # Select top 60% most variable regulons using VST method
    nhvg <- round(dim(pbmc)[1] * 0.6)
    pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = nhvg)
    hvg <- VariableFeatures(pbmc)
    print(paste0("Selected ", length(hvg), " variable regulons (60% of total)"))
    
    #####---Principal Component Analysis (PCA)----------
    # Run PCA on variable regulons
    pbmc <- RunPCA(pbmc, features = hvg, npcs = npc)
    
    #####---Replace PCA Embeddings with Differential Regulons----------
    # Instead of using actual PCA coordinates, we use the scaled expression
    # values of significant differential regulons directly as embeddings.
    # This ensures the trajectory analysis focuses on biologically meaningful
    # regulon differences between domains.
    deg.w <- as.array(t(deg.w))
    pbmc@reductions$pca@cell.embeddings <- deg.w
    print(paste0("PCA embeddings replaced with ", ncol(deg.w), " differential regulons"))
    
    #####---Non-linear Dimensionality Reduction----------
    # Use all significant differential regulons as dimensions for UMAP/tSNE
    pcs <- seq_len(length(unique(deg.s$gene)))
    
    # Run UMAP
    pbmc <- RunUMAP(pbmc, dims = pcs, n.components = 2, reduction = "pca", 
                    n.neighbors = n.neighbors)
    
    # Run tSNE
    pbmc <- RunTSNE(pbmc, reduction = "pca", dims = pcs, perplexity = perplexity)
    
    
        #####---Visualization of Dimensionality Reduction----------
    # Prepare color palettes for different annotations
    cols1 <- unique(meta.12r$colr.l)   # Lineage colors
    names(cols1) <- unique(meta.12r$lineage)
    
    # Sector colors (use metadata colors if available, otherwise generate rainbow)
    if (length(unique(meta.12r$colr.s)) == length(unique(meta.12r$sector23))) {
      cols2 <- unique(meta.12r$colr.s)
      names(cols2) <- unique(meta.12r$sector23)
    } else {
      cols2 <- rainbow(length(unique(meta.12r$sector23)))
      names(cols2) <- unique(meta.12r$sector23)
    }
    
    # Domain annotation colors
    if (length(unique(meta.12r$colr.d)) == length(unique(meta.12r$annotation))) {
      cols3 <- unique(meta.12r$colr.d)
      names(cols3) <- unique(meta.12r$annotation)
    } else {
      cols3 <- rainbow(length(unique(meta.12r$annotation)))
      names(cols3) <- unique(meta.12r$annotation)
    }
    
    # Fetch visualization data from Seurat object
    pltd <- FetchData(pbmc, vars = c("tSNE_1", "tSNE_2", "UMAP_1", "UMAP_2",
                                     "sector23", "lineage", "section", "newname",
                                     "stage", "annotation"))
    
    #####---tSNE and UMAP Plots Colored by Sector----------
    tSNE_plot <- ggplot(pltd, aes(tSNE_1, tSNE_2, color = sector23)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "tSNE")) +
      scale_colour_manual(values = cols2) +
      theme_classic()
    
    UMAP_plot <- ggplot(pltd, aes(UMAP_1, UMAP_2, color = sector23)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "UMAP")) +
      scale_colour_manual(values = cols2) +
      theme_classic()
    
    # Save combined tSNE/UMAP plot
    pdf(file = paste0(time_1, "_", time_2, "/", fn, ".PCA.tSNE.UMAP.cluster.sector.pdf"),
        width = 15, height = 8)
    p <- ggarrange(UMAP_plot, tSNE_plot, common.legend = TRUE, legend = "right")
    plot(p)
    dev.off()
    
    #####---UMAP Plots Colored by Different Annotations----------
    # Domain annotation
    DOME_plot <- ggplot(pltd, aes(UMAP_1, UMAP_2, color = annotation)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "domain")) +
      scale_colour_manual(values = cols3) +
      theme_classic()
    
    # Stage
    UMAP_stage <- ggplot(pltd, aes(UMAP_1, UMAP_2, color = stage)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "_stage")) +
      theme_classic()
    
    # Section
    UMAP_section <- ggplot(pltd, aes(UMAP_1, UMAP_2, color = section)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "_section")) +
      theme_classic()
    
    # Lineage
    UMAP_lineage <- ggplot(pltd, aes(UMAP_1, UMAP_2, color = lineage)) +
      geom_point(size = 2) +
      ggtitle(paste(fn, "_lineage")) +
      scale_colour_manual(values = cols1) +
      theme_classic()
    
    # Save combined annotation plots
    pdf(file = paste0(time_1, "_", time_2, "/", fn, ".UMAP.section.lineage.pdf"),
        width = 20, height = 12)
    p <- ggarrange(DOME_plot, UMAP_stage, UMAP_section, UMAP_lineage, 
                   common.legend = FALSE, legend = "right")
    plot(p)
    dev.off()
    
    #####---Save PCA Embeddings for KNN Analysis----------
    # Save the PCA embeddings (which contain differential regulon scores)
    # These will be used for KNN-based lineage inference
    emb <- data.frame(Embeddings(object = pbmc, reduction = "pca"))
    saveRDS(emb, file = paste0(dir, "/", time_1, "_", time_2, "_pca.rds"))
    print(paste0("Saved PCA embeddings: ", dir, "/", time_1, "_", time_2, "_pca.rds"))
    
    
  } 
  
  
    #####---KNN-based Lineage Inference----------
  # Use KNN to predict ancestor-descendant relationships between stages
  print(paste0("Performing KNN lineage inference for stage k = ", kk))
  
  time_i <- time_point[kk]
  time_j <- time_point[kk + 1]
  print(paste0("Stage transition: ", time_i, " -> ", time_j))
  
  # Add 'day' column to metadata to distinguish pre (earlier) vs nex (later) stages
  meta.12r$day <- ifelse(meta.12r$stage == time_i, "pre", "nex")
  saveRDS(meta.12r, paste0(dir, "/meta.regulon_", time_1, "_", time_2, "merged.rds"))
  
  # Prepare annotation data for KNN analysis
  anno <- meta.12r[, c("domain", "day", "stage", "stname")]
  
  # Ensure annotation matches embedding order
  if (nrow(emb) != nrow(anno)) {
    print("Error: Embedding and annotation dimensions mismatch!")
  }
  anno <- anno[rownames(emb), ]
  
  print("Domains in annotation:")
  print(unique(anno$domain))
  
  # Run KNN lineage inference with k=5 and k=10 neighbors
  # The createLineage_Knn function performs bootstrapped KNN analysis
  # to infer probabilistic connections between cell states across stages
  for (k_neigh in c(5, 10)) {
    print(paste0("Running KNN with k = ", k_neigh))
    
    # Create output directory for this k value
    dir.create(paste0("knn", k_neigh), showWarnings = FALSE)
    
    # Run lineage inference (1000 bootstrap iterations)
    res <- createLineage_Knn(emb, anno, reduction = "pca", 
                             k_neigh = k_neigh, replication_times = 1000)
    
    # Save results
    saveRDS(res, paste0(dir, "/knn", k_neigh, "/", time_i, "_", time_j, 
                        "_Knn", k_neigh, "_pca.rds"))
  }
  
}  # End of main stage loop

###############################################################################
#### Step 2: Calculate Median Edge Probabilities Across Bootstrap Iterations ####
###############################################################################

replication_times <- 1000

# Process results for both k=5 and k=10
for (k_neigh in c(5, 10)) {
  
  work_path <- paste0(dir, "/knn", k_neigh)
  res_median_umap <- list()
  
  # Process each stage pair
  for (ti in 1:(length(time_point) - 1)) {
    print(paste0("Processing stage pair: ", time_point[ti], " -> ", time_point[ti + 1]))
    
    # Load KNN results for this stage pair
    dat <- readRDS(paste0(work_path, "/", time_point[ti], "_", time_point[ti + 1],
                          "_Knn", k_neigh, "_pca.rds"))
    
    # Get state names (rows = nex states, columns = pre states)
    state_1 <- row.names(dat[[1]])  # Later stage states (nex)
    state_2 <- names(dat[[1]])       # Earlier stage states (pre)
    
    # Calculate median probability across all bootstrap iterations
    tmp_1 <- matrix(NA, nrow(dat[[1]]), ncol(dat[[1]]))
    for (i in 1:nrow(dat[[1]])) {
      for (j in 1:ncol(dat[[1]])) {
        # Collect probabilities from all iterations
        xx <- NULL
        for (k in 1:replication_times) {
          xx <- c(xx, dat[[k]][i, j])
        }
        # Compute median (ignoring NA values)
        tmp_1[i, j] <- median(xx[!is.na(xx)])
      }
    }
    
    tmp_1 <- data.frame(tmp_1)
    row.names(tmp_1) <- state_1
    names(tmp_1) <- state_2
    res_median_umap[[ti]] <- tmp_1
  }
  
  #####---Compile Edge List----------
  # Combine all stage pair results into a single edge list
  dat <- NULL
  for (i in 1:length(res_median_umap)) {
    dat <- rbind(dat, melt(as.matrix(res_median_umap[[i]])))
  }
  
  dat <- data.frame(dat)
  names(dat) <- c("nex", "pre", "prob")  # nex = later stage, pre = earlier stage, prob = probability
  
  # Parse stage and cell type from state names (format: E7.5_domain1)
  dat$pre_time <- unlist(lapply(as.vector(dat$pre), function(x) strsplit(x, "[_]")[[1]][1]))
  dat$pre_cell <- unlist(lapply(as.vector(dat$pre), function(x) strsplit(x, "[_]")[[1]][2]))
  dat$nex_time <- unlist(lapply(as.vector(dat$nex), function(x) strsplit(x, "[_]")[[1]][1]))
  dat$nex_cell <- unlist(lapply(as.vector(dat$nex), function(x) strsplit(x, "[_]")[[1]][2]))
  
  # Save edge list
  saveRDS(dat, paste0(work_path, "/edge_all_pca.regressE.rds"))
  write.csv(dat, paste0(work_path, "/edge_all_pca.regressE.csv"))
  
  #####---Log Edge Statistics----------
  # Here we use "cell state" to mean an annotated cluster at a given stage
  edgelog <- c(
    print(paste0("Total edges: ", nrow(dat))),
    print(paste0("Edges with prob > 0: ", nrow(dat[dat$prob > 0,]))),
    print(paste0("Edges with prob >= 0.2: ", nrow(dat[dat$prob >= 0.2,]))),
    print(paste0("Edges with prob >= 0.3: ", nrow(dat[dat$prob >= 0.3,]))),
    print(paste0("Edges with prob >= 0.4: ", nrow(dat[dat$prob >= 0.4,]))),
    print(paste0("Edges with prob >= 0.5: ", nrow(dat[dat$prob >= 0.5,]))),
    print(paste0("Edges with prob >= 0.7: ", nrow(dat[dat$prob >= 0.7,]))),
    print(paste0("Edges with prob >= 0.8: ", nrow(dat[dat$prob >= 0.8,]))),
    print(paste0("Total nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex)))))),
    print(paste0("Total cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))
  )
  cat(edgelog, file = paste0(work_path, "/edge.weight.log.txt"), sep = "\t", append = TRUE)
  
}

###############################################################################
#### Step 3: Edge Filtering and Sankey Plot Generation ####
###############################################################################

# Load additional trajectory data (from single-cell analysis)
mfls <- list.files(path = "../", pattern = paste0("sc.trajectory.add.csv"),
                   recursive = FALSE, full.names = TRUE)
trjadd <- read.csv(mfls, header = TRUE)
print("Additional trajectory data preview:")
head(trjadd)

# Load domain name mapping (for renaming consistency)
mfls <- list.files(path = "../", pattern = paste0("sc.domain.namechange.csv"),
                   recursive = FALSE, full.names = TRUE)
scn <- read.csv(mfls, header = FALSE)
colnames(scn) <- c("old", "new")
print("Domain name mapping preview:")
head(scn)

work_path <- dir

# Process both k=5 and k=10 results
for (k_neigh in c(5, 10)) {
  
  # Load edge list
  mfls <- list.files(path = paste0("knn", k_neigh), 
                     pattern = paste0("edge_all_pca.regressE.rds"),
                     recursive = FALSE, full.names = TRUE)
  dat <- readRDS(mfls)
  print(paste0("Loaded edge list: ", nrow(dat), " edges"))
  
  # Combine with single-cell trajectory edges
  dat <- rbind(trjadd, dat)
  print(paste0("After adding SC edges: ", nrow(dat), " edges"))
  
  # Rename domains using name mapping
  for (i in seq_len(nrow(scn))) {
    dat$nex[dat$nex == scn$old[i]] <- scn$new[i]
    dat$pre[dat$pre == scn$old[i]] <- scn$new[i]
  }
  
  #####---Filter Edges by Probability Threshold----------
  # Filter edges with probability >= 0.2
  print(paste0("Edges with prob > 0: ", nrow(dat[dat$prob > 0, ])))
  print(paste0("Edges with prob >= 0.2: ", nrow(dat[dat$prob >= 0.2, ])))
  
  top2 <- dat
  nprob <- 0.2
  topf <- top2[top2$prob >= nprob, ]
  print(paste0("Filtered edges (prob >= 0.2): ", nrow(topf)))
  
  topfp <- topf
  
  #####---Ensure All Nodes Have At Least One Edge----------
  # If some nex nodes are missing, add their highest-probability edge
  if (length(unique(topf$nex)) < length(unique(dat$nex))) {
    p2n <- data.frame()
    for (x in setdiff(unique(dat$nex), unique(topf$nex))) {
      xp <- filter(dat, nex == x)
      xpf <- xp[order(xp$prob, decreasing = TRUE)[1], ]
      p2n <- rbind(p2n, xpf)
    }
    topfpn <- rbind(topfp, p2n)
  } else {
    topfpn <- topfp
  }
  
  #####---Final Edge List Preparation----------
  x <- topfpn[, c("pre", "nex", "prob")]
  
  # Check for any missing nodes
  los.pre <- setdiff(unique(c(as.vector(dat$pre), as.vector(dat$nex))),
                     unique(c(as.vector(x$pre), as.vector(x$nex))))
  
  if (length(los.pre) > 0) {
    # Add edges for missing nodes
    xp <- filter(dat, pre == los.pre)
    xpf <- xp[order(xp$prob, decreasing = TRUE)[1], ]
    y <- rbind(topfpn, xpf)
    y <- y[, c("pre", "nex", "prob")]
    res <- y
  } else {
    res <- x
  }
  
  print(paste0("Final edge list: ", nrow(res), " edges, ", 
               length(unique(c(as.vector(res$pre), as.vector(res$nex)))), " nodes"))
  
  # Save final edge list
  write.table(res, 
              paste0(work_path, "/edge", nrow(res), "_nodes", 
                     length(unique(c(as.vector(res$pre), as.vector(res$nex)))),
                     "_knn", k_neigh, "_parent1.addSC_prob.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  #####---Generate Interactive Sankey Plot----------
  # Load domain annotation file
  d.anno <- read.csv(paste0("../domain.annotation.csv"))
  colnames(d.anno)[which(colnames(d.anno) == "domain")] <- "stg.clus"
  d.anno$domain <- paste0(d.anno$stage, "_", d.anno$annotation)
  rownames(d.anno) <- d.anno$domain
  
  # Create Sankey data frame
  snak <- data.frame(stage1 = res$pre, stage2 = res$nex, prob = res$prob)
  
  # Map edge IDs to domain indices (plotly requires 0-based indexing)
  snak$IDsource <- match(snak$stage1, d.anno$domain) - 1
  snak$IDtarget <- match(snak$stage2, d.anno$domain) - 1
  
  # Create edge labels with probability
  snak$label <- paste0(snak$stage1, ":", as.character(round(snak$prob, 3)))
  
  # Add edge colors based on source lineage
  edge_color <- data.frame(stage1 = d.anno$domain, colr.l = d.anno$colr.l)
  snak <- merge(snak, edge_color, by = 'stage1', all.x = TRUE)
  
  # Compute edge weights (exponential scaling of probability)
  snak$weight <- 100^snak$prob
  
  # Save Sankey object
  saveRDS(snak, file = paste0(work_path, "/snaky.S", length(time_point), 
                              ".regulon.edge", nrow(res), "_nodes",
                              length(unique(c(as.vector(res$pre), as.vector(res$nex)))),
                              "_knn", k_neigh, "_parent1.addSC.obj.rds"))
  
  fn <- paste0("knn", k_neigh, "_edge", nrow(res), "_nodes",
               length(unique(c(as.vector(res$pre), as.vector(res$nex)))),
               "_parent1.addSC_")
  
  #####---Create Plotly Sankey Visualization----------
  p <- plot_ly(
    type = 'sankey', orientation = 'h', arrangement = "freeform",
    
    # Node configuration
    node = list(
      label = d.anno$annotation,      # Domain annotation labels
      color = d.anno$colr.d,          # Domain colors
      pad = 50,                       # Padding between nodes
      thickness = 20,                 # Node thickness
      line = list(color = 'black', width = 0.5)  # Node border
    ),
    
    # Link configuration
    link = list(
      source = snak$IDsource,         # Source node indices
      target = snak$IDtarget,         # Target node indices
      value = snak$prob,              # Edge weights (probabilities)
      label = snak$label,             # Edge labels
      color = snak$colr.l             # Edge colors (matching source lineage)
    )
  )
  
  # Customize layout
  p <- p %>% layout(
    title = "Regulome trajectory weighted by probability",
    font = list(size = 10)
  )
  
  # Save interactive HTML plot
  htmlwidgets::saveWidget(as_widget(p), 
                          paste0(work_path, "/snakey.KNN.pca.prob.nodePosition", 
                                 nprob, "_", fn, ".html"))
  saveRDS(p, file = paste0(work_path, "/snakey.KNN.pca.prob.", nprob, "_", fn, ".obj.rds"))
  
}

# Exit R session
q()










