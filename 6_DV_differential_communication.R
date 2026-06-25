# ==============================================================================
# Signaling Pathway Differential Analysis across Dorsal-Ventral Axis
# ==============================================================================
#
# Description:
# This script analyzes CellChat signaling pathway probabilities to identify
# differential signaling between dorsal (D) and ventral (V) domains across
# multiple developmental stages. The workflow includes:
#
# Step 1: Select pathways for analysis - For each stage and source-target pair,
#         identify top pathways with balanced positive/negative D-V differences,
#         plus 7 canonical pathways (BMP, FGF, HH, RA, WNT, NODAL, HIPPO)
#
# Step 2: Extract pathway probability data from CellChat objects and 
#         calculate D-V probability differences for selected pathways
#
# Step 3: Filter for consistent pathways that show the same
#         D-V bias across all stages and sources, compute average probabilities,
#         and assign colors based on differential magnitude
#
# Step 4: Generate dot plots - Visualize signaling pathways with D (brown) and
#         V (blue) biases for Fg and Hg domains separately
#
# Input files:
#   - RDS files: CellChat objects for each stage (*_secreted_cellchat.rds)
#     located in ../E8.75/, ../E8.5/, ../E8.25/ directories
#
# Output files:
#   - PNG/PDF files: Dot plots showing D-V differential signaling for Fg domain
#     (FgDV_diff_ccc_3stage_shared_top8.*)
#   - PNG/PDF files: Dot plots showing D-V differential signaling for Hg domain
#     (HgDV_diff_ccc_3stage_shared_top8.*)
#   - CSV files: Processed data for Fg and Hg domains with probability differences
#     (FgDV_diff_ccc_3stage_shared_top8.csv, HgDV_diff_ccc_3stage_shared_top8.csv)
#
# ==============================================================================

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# ==============================================================================
# Configuration Parameters
# ==============================================================================
fg_targets <- c("FgD", "FgV")     # Foregut dorsal/ventral targets
hg_targets <- c("HgD", "HgV")     # Hindgut dorsal/ventral targets
stages <- c("E8.75", "E8.5", "E8.25")  # Developmental stages
fg_sources <- c("SM_A", "NdA")    # Foregut source cell types
hg_sources <- c("SM_P", "NdP")    # Hindgut source cell types
canonical_pw <- c("BMP", "FGF", "HH", "RA", "WNT", "NODAL", "HIPPO")

# ==============================================================================
# Helper Functions
# ==============================================================================

# Select top rows from data frame based on sorting column
# Args:
#   df: Data frame containing the data
#   value: Number of rows to select (if type="count") or percentage (if type="percent")
#   type: "percent" or "count"
#   balance: If TRUE, select balanced positive and negative values; if FALSE,
#            select based on absolute values
#   sort_column: Column name to sort by
# Returns:
#   Data frame with selected top rows
select_top_rows <- function(df, value, type = "percent", balance = FALSE, sort_column = "diff") {
  if (type != "percent" && type != "count") {
    stop("type must be 'percent' or 'count'")
  }
  
  if (!(sort_column %in% names(df))) {
    stop("sort_column must be a column in the data frame")
  }
  
  df[[sort_column]] <- as.numeric(df[[sort_column]])
  
  if (type == "percent") {
    total_count <- round(0.01 * value * nrow(df))
  } else {
    total_count <- value
  }
  
  if (balance) {
    pos_count <- sum(df[[sort_column]] > 0)
    neg_count <- sum(df[[sort_column]] < 0)
    
    half_count <- total_count / 2
    if (half_count > pos_count) half_count <- pos_count
    if (half_count > neg_count) half_count <- neg_count
    
    positive_result <- df %>%
      filter(.[[sort_column]] > 0) %>%
      arrange(desc(.[[sort_column]])) %>%
      slice(1:half_count)
    
    negative_result <- df %>%
      filter(.[[sort_column]] < 0) %>%
      arrange(.[[sort_column]]) %>%
      slice(1:half_count)
    
    result <- rbind(positive_result, negative_result)
  } else {
    result <- df %>%
      arrange(desc(abs(.[[sort_column]]))) %>%
      slice(1:total_count)
  }
  
  return(result)
}

# Generate color palette based on differential values
# Positive values (D-side stronger) get brown colors
# Negative values (V-side stronger) get blue colors
# Args:
#   values: Numeric vector of differential values
# Returns:
#   Character vector of colors corresponding to input values
colorRamp <- function(values) {
  colors <- rep(NA, length(values))
  
  d_indices <- which(values > 0)
  if (length(d_indices) > 0) {
    d_palette <- colorRampPalette(c("burlywood1", "brown3"))(length(d_indices))
    colors[d_indices[order(values[d_indices])]] <- d_palette
  }
  
  v_indices <- which(values < 0)
  if (length(v_indices) > 0) {
    v_palette <- colorRampPalette(c("#5086C4", "#B8E5FA"))(length(v_indices))
    colors[v_indices[order(values[v_indices])]] <- v_palette
  }
  
  return(colors)
}

# ==============================================================================
# Core Functions
# ==============================================================================

# Select pathways for plotting from CellChat object
# Args:
#   stage: Developmental stage (e.g., "E8.75")
#   source: Source cell type (e.g., "SM_A")
#   target: Vector of target domains (e.g., c("FgD", "FgV"))
# Returns:
#   Character vector of pathway names to include in plots
GetPlotPW <- function(stage, source, target) {
  cellchat <- readRDS(paste0("../", stage, "/", stage, "_secreted_cellchat.rds"))
  pw_prob <- data.frame(t(cellchat@netP$prob[source, target, ]))
  colnames(pw_prob) <- paste0(source, "-->", target)
  pw_prob$Diff <- pw_prob[, 1] - pw_prob[, 2]
  
  balanced_results <- select_top_rows(pw_prob, 8, type = "count", balance = TRUE, sort_column = "Diff")
  top_pw <- rownames(balanced_results)
  
  plot_pw <- unique(c(top_pw, intersect(canonical_pw, cellchat@netP$pathways)))
  return(plot_pw)
}

# Extract pathway probability data from CellChat object
# Args:
#   stage: Developmental stage
#   source: Source cell type
#   target: Vector of target domains
#   plot_pw: Vector of pathway names to extract
# Returns:
#   Data frame with Signaling, Stage, Source, D_Prob, V_Prob columns
GetMissProb <- function(stage, source, target, plot_pw) {
  cellchat <- readRDS(paste0("../", stage, "/", stage, "_secreted_cellchat.rds"))
  
  inter_pw <- intersect(plot_pw, cellchat@netP$pathways)
  pw_prob <- data.frame(t(cellchat@netP$prob[source, target, inter_pw]))
  colnames(pw_prob) <- paste0(source, "-->", target)
  pw_prob$Diff <- pw_prob[, 1] - pw_prob[, 2]
  
  # Set D probability to Diff if positive, else 0
  pw_prob[, 1] <- ifelse(pw_prob$Diff > 0, pw_prob$Diff, 0)
  # Set V probability to Diff (negative) if negative, else 0
  pw_prob[, 2] <- ifelse(pw_prob$Diff < 0, pw_prob$Diff, 0)
  
  pw_prob <- pw_prob[, -ncol(pw_prob)]  # Remove Diff column
  pw_prob$Signaling <- rownames(pw_prob)
  pw_prob$Stage <- stage
  pw_prob$Source <- source
  
  pw_prob_long <- pw_prob %>%
    mutate(D_Prob = .[[1]], V_Prob = .[[2]]) %>%
    select(Signaling, Stage, Source, D_Prob, V_Prob)
  
  return(pw_prob_long)
}

# Process domain data: filter consistent pathways and compute average probabilities
# Args:
#   df: Combined data frame from all stages and sources
#   domain: Domain name ("Fg" for foregut, "Hg" for hindgut)
# Returns:
#   Processed data frame with Signaling, Source, Prob_Diff, Target, Source.target, color
process_domain <- function(df, domain) {
  # Define source cell types for each domain
  domain_sources <- if (domain == "Fg") c("SM_A", "NdA") else c("SM_P", "NdP")
  df_domain <- df %>% filter(Source %in% domain_sources)
  
  # Check consistency: pathways should show the same D/V bias across all stages
  consistency_check <- df_domain %>%
    group_by(Signaling, Source) %>%
    summarise(
      dominant_sides = list(ifelse(D_Prob > abs(V_Prob), "D", ifelse(abs(V_Prob) > D_Prob, "V", "0"))),
      has_inconsistency = {
        sides <- unlist(dominant_sides)
        has_d <- any(sides == "D")
        has_v <- any(sides == "V")
        has_d && has_v
      },
      .groups = 'drop'
    )
  
  # Identify inconsistent pathways (show both D and V bias across stages)
  inconsistent_signaling <- consistency_check %>%
    filter(has_inconsistency) %>%
    pull(Signaling) %>%
    unique()
  
  # Keep only consistently biased pathways
  consistent_signaling <- df_domain %>%
    pull(Signaling) %>%
    unique() %>%
    setdiff(inconsistent_signaling)
  
  # Filter: pathways must be non-zero in at least 2 stages
  non_zero_check <- df_domain %>%
    filter(Signaling %in% consistent_signaling) %>%
    group_by(Signaling, Source) %>%
    summarise(
      non_zero_count = sum((D_Prob != 0) | (V_Prob != 0)),
      .groups = 'drop'
    ) %>%
    filter(non_zero_count >= 2)
  
  final_signaling <- non_zero_check %>% pull(Signaling) %>% unique()
  
  # Compute average probabilities and differences
  df_consistent <- df_domain %>%
    filter(Signaling %in% final_signaling) %>%
    group_by(Signaling, Source) %>%
    summarise(
      Mean_D_Prob = mean(D_Prob),
      Mean_V_Prob = mean(V_Prob),
      .groups = 'drop'
    ) %>%
    mutate(
      Prob_Diff = Mean_D_Prob + Mean_V_Prob,
      Target = ifelse(Prob_Diff > 0, paste0(domain, "D"), paste0(domain, "V"))
    ) %>%
    filter(Prob_Diff != 0)
  
  # Add color and source-target labels
  df_consistent$Source.target <- paste0(df_consistent$Source, "-->", df_consistent$Target)
  df_consistent$color <- colorRamp(df_consistent$Prob_Diff)
  
  # Sort Y-axis: pathways with positive Prob_Diff (D-side) first
  df_consistent <- df_consistent %>% arrange(desc(Prob_Diff), Source, Signaling)
  df_consistent$Signaling <- factor(df_consistent$Signaling, levels = unique(df_consistent$Signaling))
  
  # Sort X-axis: D targets first, then V; Nd sources first, then SM
  sorted_levels <- if (domain == "Fg") {
    c("NdA-->FgD", "SM_A-->FgD", "NdA-->FgV", "SM_A-->FgV")
  } else {
    c("NdP-->HgD", "SM_P-->HgD", "NdP-->HgV", "SM_P-->HgV")
  }
  existing_levels <- intersect(sorted_levels, unique(df_consistent$Source.target))
  df_consistent$Source.target <- factor(df_consistent$Source.target, levels = existing_levels)
  
  return(df_consistent)
}

# Generate dot plot for domain data
# Args:
#   df: Processed data frame from process_domain
#   title: Plot title
#   filename: Output file name prefix
# Returns:
#   ggplot object
plot_domain <- function(df, title, filename) {
  p <- ggplot(df, aes(x = Source.target, y = Signaling, color = color)) +
    geom_point(size = 10) +
    labs(x = NULL, title = title) +
    scale_color_identity() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_line(color = "gray"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(filename = paste0(filename, ".png"), p, width = 2.5, height = 6)
  ggsave(filename = paste0(filename, ".pdf"), p, width = 2.5, height = 6)
  
  return(p)
}

# ==============================================================================
# Main Analysis Workflow
# ==============================================================================

# Step 1: Select pathways for analysis across all stages and sources
message("=== Step 1: Selecting pathways for analysis ===")
stages_plot_pw <- c()

# Foregut sources
for (stage in stages) {
  for (fg_source in fg_sources) {
    pw <- GetPlotPW(stage, fg_source, fg_targets)
    stages_plot_pw <- c(stages_plot_pw, pw)
    message(paste("Stage", stage, "source", fg_source, ":", length(pw), "pathways"))
  }
}

# Hindgut sources
for (stage in stages) {
  for (hg_source in hg_sources) {
    pw <- GetPlotPW(stage, hg_source, hg_targets)
    stages_plot_pw <- c(stages_plot_pw, pw)
    message(paste("Stage", stage, "source", hg_source, ":", length(pw), "pathways"))
  }
}

# Get unique pathways
stages_plot_pw <- unique(stages_plot_pw)
message(paste("Total unique pathways selected:", length(stages_plot_pw)))

# Step 2: Extract pathway probability data from CellChat objects
message("\n=== Step 2: Extracting pathway probability data ===")
stages_pw_prob <- data.frame(
  Signaling = character(), 
  Stage = character(),
  Source = character(), 
  D_Prob = numeric(), 
  V_Prob = numeric()
)

# Foregut sources
for (stage in stages) {
  for (fg_source in fg_sources) {
    fg_df <- GetMissProb(stage, fg_source, fg_targets, stages_plot_pw)
    stages_pw_prob <- rbind(stages_pw_prob, fg_df)
  }
}

# Hindgut sources
for (stage in stages) {
  for (hg_source in hg_sources) {
    hg_df <- GetMissProb(stage, hg_source, hg_targets, stages_plot_pw)
    stages_pw_prob <- rbind(stages_pw_prob, hg_df)
  }
}

message(paste("Total data points extracted:", nrow(stages_pw_prob)))

# Step 3: Process domain data and generate plots
message("\n=== Step 3: Processing Foregut (Fg) domain ===")
df_fg <- process_domain(stages_pw_prob, "Fg")
message(paste("Fg consistent pathways:", nrow(df_fg)))
p_fg <- plot_domain(df_fg, "Fg (D/V)", "FgDV_diff_ccc_3stage_shared_top8")
write.csv(df_fg, "FgDV_diff_ccc_3stage_shared_top8.csv", row.names = FALSE)

message("\n=== Step 4: Processing Hindgut (Hg) domain ===")
df_hg <- process_domain(stages_pw_prob, "Hg")
message(paste("Hg consistent pathways:", nrow(df_hg)))
p_hg <- plot_domain(df_hg, "Hg (D/V)", "HgDV_diff_ccc_3stage_shared_top8")
write.csv(df_hg, "HgDV_diff_ccc_3stage_shared_top8.csv", row.names = FALSE)

# Step 4: Display plots
message("\n=== Step 5: Displaying plots ===")
print(p_fg)
print(p_hg)

# Summary
message("\n=== Analysis Complete ===")
message(paste("Fg domain:", nrow(df_fg), "pathways"))
message(paste("Hg domain:", nrow(df_hg), "pathways"))