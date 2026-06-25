# ==============================================================================
# Spatial Analysis of Secretory Gene Signaling and Cell-Cell Communication
# ==============================================================================
#
# Description:
# This script performs spatial analysis of secretory gene expression in mouse
# embryo E7.5 endoderm. The workflow includes:
#   1. Loading 2D spatial coordinates with gene expression scores
#   2. Binarizing signaling scores using Otsu's method
#   3. Identifying high-signaling cells and filtering sparse spatial points
#   4. Finding neighboring cells around high-signaling cells
#   5. Performing CellChat analysis for cell-cell communication inference
#   6. Visualizations
#
# Input files:
#   - CSV file: 2D coordinates with expression scores (*_ST_2D_coord_cbind_mapped_SC_cell_type_meta_score.csv)
#   - CSV file: Whole embryo tangram-mapped coordinates (E7.5_ST_coord_cbind_mapped_SC_cell_type_meta.csv)
#   - CSV file: Single cell count data (E7.5_jhon_sc_whole_embryo_count.csv)
#   - CSV file: Color sheet for cell type visualization (extend_altas_colorsheet.csv)
#
# Output files:
#   - CSV file: Selected high-signaling cells (E7.5_Endo_signaling_high_cells_selected.csv)
#   - CSV file: Neighbor mapping for high-signaling cells (*_to_neighbor_*.csv)
#   - CSV file: Selected neighbor types (*_neighbor_type_select.csv)
#   - CSV file: Neighbor cell type distribution (*_neighbor_celltype_*.csv)
#   - CSV file: Filtered CellChat pathway results (*_CCC_filtered_netp.csv)
#   - PDF file: Spatial plot of selected high-signaling cells (*_signaling_high_cells_selected.pdf)
#   - PDF file: Neighbor cell type bar plot (*_neighbor_celltype_bar_*.pdf)
#   - PDF file: CellChat LR communication plot (*_CCC_LR.pdf)
#   - PDF file: CellChat pathway communication plot (*_CCC_pathway.pdf)
#   - RDS file: CellChat object (*_cellchat.rds)
#
# ==============================================================================

# Load required packages
library(Seurat)
library(AUCell)
library(ggplot2)
library(dplyr)
library(viridis)
library(dbscan)
library(CellChat)
library(patchwork)
library(reshape2)

# Import theme_academic from help_code
source("help_code/Plot.Functions.R")

# Load color sheet
colorsheet <- read.csv("../extend_altas_colorsheet.csv")
my_color <- colorsheet$color
names(my_color) <- colorsheet$cell_type

# Load 2D spatial coordinates with expression scores
coord_2d <- read.csv("/home/cui_guizhong/cgz/LCM/organogenesis/3Dsampling/3dto2d/E7.5/E7.5_Endoderm_ST_2D_coord_cbind_mapped_SC_cell_type_meta_score.csv", row.names = 1)

# ==============================================================================
# Otsu Binarization Function
# ==============================================================================
# Implements Otsu's method for automatic threshold selection based on maximizing
# between-class variance in histogram
# Args:
#   x: Numeric vector to binarize
#   na.rm: Logical, whether to remove NA values (default: TRUE)
# Returns:
#   Binary vector (0/1) with attributes 'threshold' and 'method'
# ==============================================================================
otsu_binarize <- function(x, na.rm = TRUE) {
  x_clean <- x[!is.na(x)]
  
  if (all(x_clean %in% c(0, 1))) {
    result <- as.numeric(x)
    attr(result, "threshold") <- 0.5
    attr(result, "method") <- "otsu"
    return(result)
  }
  
  hist_info <- hist(x_clean, breaks = 100, plot = FALSE)
  total <- sum(hist_info$counts)
  sum_total <- sum(hist_info$counts * hist_info$mids)
  cumsum_counts <- cumsum(hist_info$counts)
  cumsum_values <- cumsum(hist_info$counts * hist_info$mids)
  
  between_var <- (cumsum_values / cumsum_counts - 
                  (sum_total - cumsum_values) / (total - cumsum_counts))^2 * 
                 cumsum_counts * (total - cumsum_counts) / total^2
  
  thresh <- hist_info$mids[which.max(between_var)]
  result <- as.numeric(x > thresh)
  result[is.na(x)] <- NA
  
  attr(result, "threshold") <- thresh
  attr(result, "method") <- "otsu"
  
  return(result)
}

# Apply Otsu binarization to secret_ligand scores
coord_2d$bin_otsu <- otsu_binarize(coord_2d$secret_ligand)

# ==============================================================================
# Step 1: Visualize secret_ligand distribution and binarized results
# ==============================================================================
p1 <- ggplot(coord_2d, aes(x = x_larp, y = y_larp, color = secret_ligand)) +
  geom_point(size = 0.8) +
  scale_color_viridis(option = "D") +
  ggtitle("Secret ligand") +
  coord_fixed(ratio = 1) +
  theme_academic()

p2 <- ggplot(coord_2d, aes(x = x_larp, y = y_larp, color = factor(bin_otsu))) +
  geom_point(size = 0.8) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  ggtitle("Otsu binarized") +
  coord_fixed(ratio = 1) +
  theme_academic()

print(p1 + p2 + plot_layout(ncol = 1))

# ==============================================================================
# Step 2: Identify high-signaling cells (bin_otsu == 1) 
# ==============================================================================
keep_celltypes <- coord_2d %>%
  filter(bin_otsu == 1) %>%
  count(cell_type, name = "count") %>%
  mutate(percent = count / sum(count) * 100) %>%
  filter(percent > 10) %>%
  pull(cell_type)

plot_filtered <- coord_2d %>%
  filter(bin_otsu == 1, cell_type %in% keep_celltypes)

# ==============================================================================
# Step 3: Spatial density filtering using KNN to remove sparse points
# ==============================================================================
coords <- as.matrix(plot_filtered[, c("x_larp", "y_larp")])
knn_dist <- dbscan::kNN(coords, k = 30)
avg_dist <- rowMeans(knn_dist$dist)
threshold <- quantile(avg_dist, 0.85)
plot_filtered_dense <- plot_filtered[avg_dist < threshold, ]

# Plot filtered cells
p_source_filtered <- ggplot(plot_filtered_dense, aes(x = x_larp, y = y_larp, color = cell_type)) +
  geom_point(size = 3) +
  scale_color_manual(values = my_color) +
  coord_fixed(ratio = 1) +
  theme_academic() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom")

print(p_source_filtered)
ggsave("E7.5_Endo_signaling_high_cells_selected.pdf", p_source_filtered, width = 15, height = 6)
write.csv(plot_filtered_dense, "E7.5_Endo_signaling_high_cells_selected.csv")

# ==============================================================================
# Step 4: Find neighboring cells around high-signaling cells
# ==============================================================================
whole_embryo_tangram_coord <- read.csv("E7.5_ST_coord_cbind_mapped_SC_cell_type_meta.csv", row.names = 1)

# Calculate global average nearest neighbor distance in 3D
coords_3d_all <- as.matrix(whole_embryo_tangram_coord[, c("x", "y", "z")])
nn_dist_all <- dbscan::kNNdist(coords_3d_all, k = 1)
avg_global_dist <- mean(nn_dist_all)

# Search radius = 4 * average nearest neighbor distance
dist_factor <- 4
search_radius <- avg_global_dist * dist_factor

# Prepare coordinates for neighbor search
target_coords <- as.matrix(plot_filtered_dense[, c("x", "y", "z")])
ref_coords <- coords_3d_all

# Find neighbors within search radius
nn_res <- dbscan::frNN(ref_coords, eps = search_radius, query = target_coords)
neighbor_indices <- nn_res$id

# Build neighborhood mapping table
target_ids <- plot_filtered_dense$unique_id
target_types <- plot_filtered_dense$cell_type

cell_type_stats <- lapply(seq_along(neighbor_indices), function(i) {
  idx <- neighbor_indices[[i]]
  current_id <- target_ids[i]
  current_type <- target_types[i]
  
  if (length(idx) == 0) return(NULL)
  
  neighbors_meta <- whole_embryo_tangram_coord[idx, ]
  neighbors_meta <- neighbors_meta[neighbors_meta$unique_id != current_id, ]
  
  if (nrow(neighbors_meta) == 0) return(NULL)
  
  data.frame(
    central_id = current_id,
    central_type = current_type,
    neighbor_id = neighbors_meta$unique_id,
    neighbor_type = neighbors_meta$cell_type,
    stringsAsFactors = FALSE
  )
})

neighborhood_df <- bind_rows(cell_type_stats)

# Save neighborhood mapping
outfile_prefix <- "E7.5_Endo_signaling_high_cells_selected_to_neighbor"
outfile_tag <- paste0("dist_factor_", dist_factor)
write.csv(neighborhood_df, paste0(outfile_prefix, "_", outfile_tag, ".csv"), row.names = FALSE)

# ==============================================================================
# Step 5: Analyze neighbor cell type composition
# ==============================================================================
# Global cell type distribution among neighbors
global_cell_type_dist <- neighborhood_df %>%
  group_by(neighbor_type) %>%
  summarise(total_count = n()) %>%
  mutate(percentage = (total_count / sum(total_count)) * 100) %>%
  arrange(desc(percentage))

# Plot neighbor cell type composition
ratio_plot <- ggplot(global_cell_type_dist, aes(x = reorder(neighbor_type, -percentage), y = percentage, fill = neighbor_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_color) +
  theme_academic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Composition of Neighborhood Cell Types",
    subtitle = paste0("Total active points analyzed: ", length(unique(neighborhood_df$central_id))),
    x = "Neighbor Cell Type",
    y = "Percentage (%)"
  )

print(ratio_plot)
ggsave(paste0(outfile_prefix, "_neighbor_celltype_bar_", outfile_tag, ".pdf"), ratio_plot, width = 10, height = 6)
write.csv(global_cell_type_dist, paste0(outfile_prefix, "_neighbor_celltype_", outfile_tag, ".csv"), row.names = FALSE)

# ==============================================================================
# Step 6: Select top cell types contributing 80% of neighbors
# ==============================================================================
target_categories <- global_cell_type_dist %>%
  arrange(desc(percentage)) %>%
  mutate(cum_pct = cumsum(percentage)) %>%
  filter(lag(cum_pct, default = 0) < 80) %>%
  pull(neighbor_type)

neighbor_type_selected_neighborhood_df <- neighborhood_df %>%
  filter(neighbor_type %in% target_categories)

write.csv(neighbor_type_selected_neighborhood_df,
          "E7.5_Endo_signaling_high_cells_selected_to_neighbor_dist_factor_4_neighbor_type_select.csv",
          row.names = FALSE)

# ==============================================================================
# Step 7: CellChat Analysis
# ==============================================================================
# Load count data and build meta mapping
count_data <- read.csv("E7.5_jhon_sc_whole_embryo_count.csv", row.names = 1)

# Extract source and target cell IDs
final_target_unique_ids <- neighborhood_df %>%
  filter(neighbor_type %in% target_categories) %>%
  pull(neighbor_id) %>%
  unique()

final_source_unique_ids <- neighborhood_df %>%
  pull(central_id) %>%
  unique()

# Build meta mapping from tangram data
tangram_data <- read.csv("E7.5_whole_embryo_tangram_mapped_sp_cell_one2one.csv")
meta_mapping <- tangram_data %>%
  select(SC, unique_id) %>%
  distinct()

# Get SC names for source and target cells
target_sc_names <- meta_mapping %>%
  filter(unique_id %in% final_target_unique_ids) %>%
  pull(SC) %>%
  unique()

source_sc_names <- meta_mapping %>%
  filter(unique_id %in% final_source_unique_ids) %>%
  pull(SC) %>%
  unique()

# Prepare count matrix and metadata for CellChat
all_sc_to_use <- unique(c(target_sc_names, source_sc_names))
available_sc <- intersect(all_sc_to_use, colnames(count_data))
count_sub <- as.matrix(count_data[, available_sc])

# Build CellChat metadata
meta_cc <- data.frame(
  SC_ID = available_sc,
  row.names = available_sc
)

# Map cell_type to SC
sc_celltype_map <- meta_mapping %>%
  left_join(whole_embryo_tangram_coord[, c("unique_id", "cell_type")], by = "unique_id") %>%
  filter(SC %in% available_sc) %>%
  distinct(SC, .keep_all = TRUE) %>%
  select(SC, cell_type)

meta_cc$cell_type <- sc_celltype_map$cell_type[match(rownames(meta_cc), sc_celltype_map$SC)]

# Label cells as Source, Target, Both, or Other
meta_cc$group <- case_when(
  rownames(meta_cc) %in% source_sc_names & rownames(meta_cc) %in% target_sc_names ~ "Both",
  rownames(meta_cc) %in% source_sc_names ~ "Source",
  rownames(meta_cc) %in% target_sc_names ~ "Target",
  TRUE ~ "Other"
)

# Create CellChat object
cellchat <- createCellChat(object = count_sub, meta = meta_cc, group.by = "cell_type")

# Load CellChat database for mouse and filter for secreted signaling
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use

# Run CellChat pipeline
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 2)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "E7.5_Endo_signaling_high_cells_selected_to_neighbor_cellchat.rds")

# ==============================================================================
# Step 8: Extract and visualize CellChat results
# ==============================================================================
# Identify source and target cell types
src_types <- unique(meta_cc$cell_type[meta_cc$group %in% c("Source", "Both")])
tgt_types <- unique(meta_cc$cell_type[meta_cc$group %in% c("Target", "Both")])

# Extract ligand-receptor communication network
df_net <- subsetCommunication(cellchat,
                              sources.use = src_types,
                              targets.use = tgt_types) %>%
  filter(source != target)

df_net <- df_net %>% arrange(desc(prob))
df_net$source.target <- paste0(df_net$source, "-->", df_net$target)

# Assign p-value categories for visualization
df_net$pval <- with(df_net,
  ifelse(pval > 0.05, 1,
  ifelse(pval > 0.01, 2, 3)))

# Plot LR communication network
font.size <- 14
colors <- c("#FFD700", "#FF0000")

g_lr <- ggplot(df_net, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
  geom_point(pch = 16) +
  theme_linedraw() +
  theme(panel.grid.major = element_line(color = "grey90")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom") +
  scale_radius(range = c(min(df_net$pval) * 1.5, max(df_net$pval) * 1.5),
               breaks = sort(unique(df_net$pval)),
               labels = c("p > 0.05", "0.01 < p ≤ 0.05", "p ≤ 0.01"),
               name = "p-value") +
  scale_colour_gradientn(colors = colors,
                        limits = c(quantile(df_net$prob, 0, na.rm = TRUE), quantile(df_net$prob, 1, na.rm = TRUE)),
                        breaks = c(quantile(df_net$prob, 0, na.rm = TRUE), quantile(df_net$prob, 1, na.rm = TRUE)),
                        labels = c("min", "max")) +
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
  theme(text = element_text(size = font.size),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6))

print(g_lr)
ggsave("E7.5_Endo_signaling_high_cells_selected_to_neighbor_CCC_LR.pdf", g_lr)

# Extract pathway-level communication network
netP <- cellchat@netP
prob <- netP$prob
netP <- reshape2::melt(prob, value.name = "prob")
colnames(netP)[1:3] <- c("source", "target", "pathway")

filtered_netp <- netP %>%
  filter(prob > 0,
         source %in% src_types,
         target %in% tgt_types,
         source != target)

write.csv(filtered_netp, "E7.5_Endo_signaling_high_cells_selected_to_neighbor_CCC_filtered_netp.csv", row.names = FALSE)

filtered_netp$source.target <- paste0(filtered_netp$source, "-->", filtered_netp$target)
filtered_netp$pathway <- factor(filtered_netp$pathway, levels = unique(filtered_netp$pathway))

# Plot pathway communication network
g_pathway <- ggplot(filtered_netp, aes(x = source.target, y = pathway, color = prob)) +
  geom_point(pch = 16, size = 6) +
  theme_linedraw() +
  theme(panel.grid.major = element_line(color = "grey90")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "bottom") +
  scale_colour_gradientn(colors = colors,
                        na.value = NA,
                        limits = c(quantile(filtered_netp$prob, 0, na.rm = TRUE), quantile(filtered_netp$prob, 1, na.rm = TRUE)),
                        breaks = c(quantile(filtered_netp$prob, 0, na.rm = TRUE), quantile(filtered_netp$prob, 1, na.rm = TRUE)),
                        labels = c("min", "max")) +
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) +
  theme(text = element_text(size = font.size),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        plot.margin = margin(l = 20, unit = "mm"))

print(g_pathway)
ggsave("E7.5_Endo_signaling_high_cells_selected_to_neighbor_CCC_pathway.pdf", g_pathway)