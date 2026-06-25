# ==============================================================================
# Regulon-based Domain Clustering and Module Analysis
# ==============================================================================

# Description:
# This script performs domain-level regulon analysis through two main steps:
#   1. Perform hierarchical clustering and visualize domain relationships as a tree
#   2. Analyze module-trait correlations between regulon modules and domain groups
#
# Input files:
#   - AUC_ave_by_domain_E5.5_E8.75.csv: Averaged AUC scores by domain
#   - metadata_E3.5_E8.75.csv: Sample metadata for domain colors
#   - domain_groups_E5.5_E8.75.csv: Domain grouping information
#   - regulon_modules.csv: Regulon module assignments
#
# Output files:
#   - domain_hclust_tree.pdf: Circular tree visualization of domains
#   - reg_module_domain_group_relationship.pdf: Module-trait correlation heatmap
#
# ==============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggtree)

# Set random seed for reproducibility
set.seed(124)

# ==============================================================================
# 1: Hierarchical Clustering and Tree Visualization
# ==============================================================================

# Read domain-averaged AUC scores
score_by_domain <- read.csv("../data/AUC_ave_by_domain_E5.5_E8.75.csv",
                            row.names = 1,
                            check.names = FALSE)

# Read metadata for domain annotations
meta <- read.csv("../data/metadata_E3.5_E8.75.csv",
                 header = TRUE,
                 row.names = 1)

# Create unique domain annotation dataframe
df_unique_domain <- meta %>%
  select(stage, annotation, colr.d) %>%
  distinct() %>%
  mutate(domain = annotation) %>%
  select(-annotation)

# Set row names as stage_domain format
rownames(df_unique_domain) <- paste0(df_unique_domain$stage, "_", df_unique_domain$domain)

# Match annotations to score_by_domain columns
annotation_col <- df_unique_domain[colnames(score_by_domain), ]

# Convert to factors
annotation_col$stage <- factor(annotation_col$stage)
annotation_col$domain <- factor(annotation_col$domain)

# Perform hierarchical clustering with pheatmap
out <- pheatmap::pheatmap(score_by_domain,
                          scale = "row",
                          clustering_method = "ward.D",
                          silent = TRUE)

# Extract column tree
hc <- out$tree_col

# Create circular tree plot with ggtree
p <- ggtree(hc, layout = "circular", branch.length = "full")
d <- data.frame(label = hc$labels, 
                colr.d = annotation_col[hc$labels, "colr.d"])
tree <- p %<+% d + 
  geom_tippoint(aes(color = colr.d), size = 3) +
  scale_color_identity() +
  geom_tiplab(aes(label = label, color = colr.d), size = 4, offset = 10)

# Save tree plot
ggsave(tree, filename = "../figure/domain_hclust_tree.pdf", 
       width = 11, height = 10)

# ==============================================================================
# 2: Module-Trait Correlation Analysis
# ==============================================================================

# Read clustering annotation metadata
annotation_col_meta <- read.csv("../data/domain_groups_E5.5_E8.75.csv",
                                row.names = 1)

# Create trait data (group indicators)
df_wide <- annotation_col_meta %>%
  group_by(group) %>%
  mutate(indicator = 1) %>%
  ungroup()

datTrait <- df_wide %>%
  pivot_wider(names_from = group, values_from = indicator, values_fill = 0)
datTrait <- datTrait[, -1]  # Remove group column

# Read module assignments
module_df <- read.csv("../data/regulon_modules.csv",
                      row.names = 1,
                      check.names = FALSE)

# Prepare input matrix (domains as rows, modules as columns)
input_mat <- t(score_by_domain[rownames(module_df), rownames(annotation_col_meta)])

# Calculate module eigengenes using WGCNA
nSamples <- nrow(input_mat)
me_list <- moduleEigengenes(input_mat, module_df$csi_module)
MEs0 <- me_list$eigengenes
MEs <- orderMEs(MEs0)

# Calculate module-trait correlations
moduleTraitCor <- WGCNA::cor(MEs, datTrait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Get group and module levels
group_level <- unique(annotation_col_meta$group)
module_level <- paste0("ME", unique(module_df$csi_module))

# Reorder correlation matrix
reorder_moduleTraitCor <- moduleTraitCor[module_level, group_level]
reorder_moduleTraitPvalue <- moduleTraitPvalue[module_level, group_level]

# Convert p-values to significance symbols
convert_to_symbols <- function(x) {
  if (x > 0.05) {
    return("")
  } else if (x > 0.01 && x <= 0.05) {
    return("*")
  } else if (x > 0.001 && x <= 0.01) {
    return("**")
  } else if (x <= 0.001) {
    return("***")
  }
}

# Apply symbol conversion to p-value matrix
symbol_reorder_moduleTraitPvalue <- apply(reorder_moduleTraitPvalue, c(1, 2), convert_to_symbols)

# Define color vector for heatmap
color_vector <- c(colorRampPalette(c("darkgreen", "white"))(100), 
                  colorRampPalette(c("white", "red"))(100))

# Create color function for heatmap
col_fun <- colorRamp2(seq(-1, 1, length.out = length(color_vector)), color_vector)

# Create heatmap with significance symbols
ht <- Heatmap(
  reorder_moduleTraitCor,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  width = unit(ncol(reorder_moduleTraitCor) * 8, "mm"),
  height = unit(nrow(reorder_moduleTraitCor) * 8, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    symbol <- symbol_reorder_moduleTraitPvalue[i, j]
    if (symbol != "") {
      grid.text(
        label = symbol,
        x = x, y = y,
        gp = gpar(fontsize = 14, fontface = "bold", col = "black")
      )
    }
  }
)

# Save heatmap as PDF
pdf("../figure/reg_module_domain_group_relationship.pdf", 
    width = 8, height = 6)
draw(ht)
dev.off()

