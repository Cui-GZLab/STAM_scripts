## ==========================================================================================
## STAM Plot Functions Library
## ==========================================================================================
## This script contains a collection of R utility and visualization functions used for
## analyzing and visualizing organogenesis data from LCM experiments. Key functions include:
## 
## 1. Utility functions: getmode, label_position
## 2. Data processing: conditionAve, dis2fram, multmerge
## 3. Visualization: corn.p (section plot), ggscatter, ggdehist, p.heat (heatmap), p.venn
## 4. Lineage inference: createLineage_Knn (KNN-based trajectory reconstruction)
## 
## Author: guizhong cui
## Date: Fri Sep 16 14:59:54 2022
## ==========================================================================================


#### Utility Functions for Clustree Plot ####

# Calculate the mode (most frequent value) of a vector
# Args:
#   v: input vector
# Returns:
#   The most frequent value in the vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Determine the label position for clustree visualization
# If all labels are the same, return that label; otherwise return "m." followed by the mode
# Args:
#   labels: vector of labels
# Returns:
#   Position label string
label_position <- function(labels) {
  if (length(unique(labels)) == 1) {
    position <- as.character(unique(labels))
  } else {
    position <- paste0("m.",getmode(labels))
  }
  return(position)
}



####--Corn Plot: Spatial Visualization of LCM Sections-------------

# Generate a spatial plot showing the distribution of LCM-captured tissue sections
# This plot visualizes the spatial arrangement of samples along the proximal-distal axis
# Args:
#   stage: Developmental stage label for the plot caption
#   phomd: Metadata data frame containing spatial coordinates (X, Y), section info, 
#          sector annotations (sector23), lineage (lineage), and sample type (type)
#          Required columns: X, Y, section, sector23, lineage, type, colr.s, colr.l
# Returns:
#   ggplot object showing the spatial arrangement of LCM sections
corn.p <- function(stage,phomd){
  pcpho <- phomd
  # which(is.na(pcpho))
  
  # Extract color mappings for lineages
  cols1 <-unique( pcpho$colr.l)
  names(cols1) <- unique(pcpho$lineage)
  cols1 
  
  # Extract color mappings for sectors (sector23)
  cols2 <- unique(pcpho$colr.s)
  names(cols2) <- unique(pcpho$sector23)
  cols2
  
  # Extract shape mappings for sample types
  shp <-as.integer(unique(pcpho$type) ) 
  names(shp) <- as.character(unique(pcpho$type))
  shp
  
  # Generate Y-axis labels for sections (formatted as "S1", "S3", etc.)
  ysec <- paste0("S",sort(unique(pcpho$section))*2-1 ) 
  names(ysec) <- sort(unique(pcpho$Y)) 
  ysec
  
  # Separate samples by type: type 16/20 (hollow shapes) vs type 21/23 (filled shapes)
  pd1 <- pcpho[pcpho$type %in% c("16","20") ,]
  pd2 <- pcpho[pcpho$type %in% c("21","23"),]
  
  # Create ggplot with spatial points
  pp<-ggplot() +
    # Plot hollow-shaped samples (type 16/20)
    geom_point(data = pd1,aes(X,Y,shape=type,color=sector23), size=6) +
    # Plot filled-shaped samples (type 21/23)
    geom_point(data = pd2,aes(X,Y,shape=type,fill=sector23), size=6) + 
    labs(caption=stage)+
    xlab(paste0("Distoal")) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA),
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          plot.caption = element_text(hjust=0.5, size=rel(3)),
          legend.position = "none")  +
    scale_colour_manual(values =cols2 ) +
    scale_shape_manual(values=shp)+
    scale_fill_manual(values =cols2 ) +
    scale_y_continuous(name = "Section", breaks = seq_len(length(ysec)), labels = ysec)
  return(pp)
}


####----GGally Custom Plot Functions for ggpairs--------

# Custom scatter plot function for GGally::ggpairs
# Displays points with density-based coloring (darker = higher density)
# Args:
#   data: Data frame
#   mapping: ggplot aesthetic mapping
# Returns:
#   ggplot scatter plot with point density coloring
ggscatter <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  df <- data.frame(x = x, y = y)
  sp1 <- ggplot(df, aes(x=x, y=y)) +
    geom_point() +
    geom_pointdensity(adjust = 1) +  # Color points by local density
    scale_color_viridis()  # Viridis color scale for density
  # geom_abline(intercept = 0, slope = 1, col = 'darkred')
  return(sp1)
}

# Custom density histogram function for GGally::ggpairs
# Displays histogram with overlaid density curve
# Args:
#   data: Data frame
#   mapping: ggplot aesthetic mapping
# Returns:
#   ggplot histogram with density curve
ggdehist <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  df <- data.frame(x = x)
  dh1 <- ggplot(df, aes(x=x)) +
    geom_histogram(aes(y=..density..), bins = 5, fill = 'lightblue', color='darkred', alpha=.4) +
    geom_density(aes(y=..density..)) +  # Overlay density curve
    theme_minimal()
  return(dh1)
}


#####---KNN-based Lineage Inference Algorithm (adapted from TOME)---
#####################################################
### Function: Finding Ancestor-Decendant Relationships ###
#####################################################

# Reconstruct cell/tissue lineage relationships between two consecutive stages using KNN
# This function performs bootstrapped KNN analysis to infer probabilistic connections
# between cell states in the earlier stage ("pre") and later stage ("nex")
#
# Args:
#   emb: Matrix of cell embeddings (rows = cells, columns = dimensions from PCA/UMAP)
#   pd: Metadata data frame with required columns: 'domain' (cell state annotation), 
#       'day' (stage indicator: "pre" for earlier, "nex" for later)
#   reduction: Name of dimensionality reduction method (default: "pca")
#   replication_times: Number of bootstrap iterations (default: 500)
#   removing_cells_ratio: Proportion of cells to randomly remove in each iteration (default: 0.2)
#   k_neigh: Number of nearest neighbors to consider (default: 5)
#
# Returns:
#   List of matrices (length = replication_times), each matrix contains:
#     - Rows: Cell states in "nex" stage
#     - Columns: Cell states in "pre" stage  
#     - Values: Probability of connection (proportion of KNN hits)
#
# Algorithm:
# 1. For each bootstrap iteration:
#    a. Randomly subsample cells (removing removing_cells_ratio of cells)
#    b. Separate cells into "pre" (earlier stage) and "nex" (later stage) groups
#    c. For each cell in "nex", find k nearest neighbors in "pre"
#    d. Count how many neighbors belong to each "pre" cell state
#    e. Normalize to get probability distribution over "pre" states
# 2. Return list of probability matrices across all iterations
createLineage_Knn <- function(emb, pd, reduction="pca", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5){
  
  print(dim(emb))
  
  # Validate required metadata columns
  if(!"domain" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no domain or day in pd")}
  
  # Validate that cell names match between embeddings and metadata
  if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
  
  # Use 'domain' column as cell state identifier
  pd$state = pd$domain
  
  # Initialize list to store results from each bootstrap iteration
  res = list()
  
  rep_i = 1
  
  # Perform bootstrap iterations
  while(rep_i < (replication_times+1)){
    
    # Step 1: Randomly subsample cells (keep 80% by default)
    sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))
    emb_sub = emb[sampling_index,]
    pd_sub = pd[sampling_index,]
    
    # Step 2: Separate cells into pre (earlier) and nex (later) stage groups
    irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]  # Earlier stage embeddings
    irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]  # Later stage embeddings
    pd_sub1 <- pd_sub[pd_sub$day == "pre",]  # Earlier stage metadata
    pd_sub2 <- pd_sub[pd_sub$day == "nex",]  # Later stage metadata
    
    # Step 3: Adjust k_neigh if some pre states have too few cells
    pre_state_min = min(table(as.vector(pd_sub1$state)))
    
    # Reduce k_neigh if minimum cells per state is less than k_neigh but >= 3
    if (pre_state_min < k_neigh & pre_state_min >= 3){
      k_neigh = pre_state_min
      print(k_neigh)
    }
    
    # Skip iteration if any pre state has fewer than 3 cells (insufficient data)
    if (pre_state_min < 3){
      next
    }
    
    # Step 4: Find k nearest neighbors for each nex cell in pre cell space
    neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
    
    # Step 5: Map neighbor indices to cell states
    tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
    for(i in 1:k_neigh){
      tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
    }
    
    # Get unique state names for pre and nex stages
    state1 <- names(table(as.vector(pd_sub1$state)))  # Pre stage states
    state2 <- names(table(as.vector(pd_sub2$state)))  # Nex stage states
    
    # Step 6: Count neighbor state occurrences and compute probabilities
    tmp2 <- matrix(NA,length(state2),length(state1))
    for(i in 1:length(state2)){
      # Collect all neighbor states for cells in current nex state
      x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
      for(j in 1:length(state1)){
        # Count how many times each pre state appears as neighbor
        tmp2[i,j] <- sum(x==state1[j])
      }
    }
    
    # Normalize rows to get probability distribution
    tmp2 <- tmp2/apply(tmp2,1,sum)
    tmp2 <- data.frame(tmp2)
    row.names(tmp2) = state2  # Rows: nex states
    names(tmp2) = state1      # Columns: pre states
    
    # Store result for this iteration
    res[[rep_i]] = tmp2
    
    rep_i = rep_i + 1
    
  }
  
  return(res)
}

# --- Academic Theme for Publication-Quality Plots ---
theme_academic <- function(base_size = 10, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.5),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text  = element_text(colour = "black", size = rel(0.9)),
      axis.title = element_text(colour = "black", size = rel(1.1)),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(size = rel(0.8)),
      legend.title      = element_text(size = rel(0.9), face = "bold"),
      legend.position   = "right",
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(0.9), face = "bold", margin = margin(b = 5)),
      plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.2), margin = margin(b = 10))
    )
}
