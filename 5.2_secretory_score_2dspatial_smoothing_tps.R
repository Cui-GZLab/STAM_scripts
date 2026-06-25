# ==============================================================================
# Spatial Smoothing with Thin Plate Spline (TPS)
# ==============================================================================
#
# Description:
# This script performs spatial smoothing on 2D projected coordinates using
# Thin Plate Spline (TPS) interpolation. It processes three germ layers 
# (Endoderm, Mesoderm, Ectoderm) separately, integrating 2D spatial coordinates
# with gene expression scores. The workflow includes:
#   1. Nearest neighbor distance calculation to determine spatial resolution
#   2. Data integration of coordinates and expression scores
#   3. TPS interpolation with GCV auto-selected lambda or manual lambda
#   4. Distance-based filtering to remove interpolated values far from original points
#   5. Multi-panel visualization generation
#
# Input files:
#   - CSV files: 2D projected coordinates for each germ layer (x_larp, y_larp columns)
#   - CSV file: 3D_score file containing expression scores (SC, cell_type, secret_ligand)
#   - CSV file: Color sheet for cell_type visualization
#
# Output files:
#   - CSV files: Integrated data for each germ layer (*_ST_2D_coord_cbind_mapped_SC_cell_type_meta_score.csv)
#   - PDF files: Visualization for each target column (*_2D_visualizations_tps_*.pdf)
#   - PDF files: Summary of filtered plots for each germ layer (*_all_score_2D_filtered_summary_tps_*.pdf)
#
# ==============================================================================


#library(spatstat)
library(ggplot2)
library(fields)
library(viridis)

######################################################################
# Method 1: Basic implementation (no external packages)
# Calculates average nearest neighbor distance for spatial resolution estimation
calc_avg_nnd <- function(data) {
  # 1. Extract x, y coordinates and convert to matrix (for distance calculation)
  coords <- as.matrix(data[, c("x", "y")])
  n <- nrow(coords)  # Total number of points

  # Edge case: Fewer than 2 points, return NA (no nearest neighbor)
  if (n < 2) {
    warning("Number of points must be >= 2 to calculate nearest neighbor distance!")
    return(NA)
  }

  # 2. Calculate Euclidean distance matrix for all point pairs
  dist_matrix <- dist(coords, method = "euclidean")

  # 3. Convert to matrix format (dist returns lower triangular vector by default)
  dist_matrix <- as.matrix(dist_matrix)

  # 4. Find nearest neighbor distance for each point (exclude self: set diagonal to Inf, then take min of each row)
  diag(dist_matrix) <- Inf  # Set self-distance to infinity (not considered for minimization)
  nearest_dist <- apply(dist_matrix, 1, min)  # Min value per row = nearest neighbor distance for that point

  # 5. Calculate average nearest neighbor distance
  avg_nnd <- mean(nearest_dist)
  return(list(
    near.dist = nearest_dist,
    ave.dist = avg_nnd
  ))
}

# Data integration function
# Merges 2D coordinates with expression scores by unique sample identifier
integrate_data <- function(coord_file, score_file, stage) {
  # Read 2D coordinate file
  coord <- read.csv(coord_file)
  message("Reading 2D coordinate file: ", coord_file)
  message("2D coordinate data rows: ", nrow(coord))
  
  # Add unique_id column
  coord$unique_id <- paste(stage, coord$name, coord$ID_unit, sep = "_")
  rownames(coord) <- coord$unique_id
  message("Unique ID column added")
  
  # Read 3D_score file
  score <- read.csv(score_file)
  message("Reading 3D_score file: ", score_file)
  message("3D_score data rows: ", nrow(score))
  
  # Ensure score file also has unique_id column
  if (!"unique_id" %in% colnames(score)) {
    score$unique_id <- paste(stage, score$name, score$ID_unit, sep = "_")
  }
  
  # Keep only required columns: SC, cell_type, secret_ligand
  score_subset <- score[, c("unique_id", setdiff(colnames(score), colnames(coord)))]
  
  # Integrate data by unique_id
  integrated <- merge(coord, score_subset, by = "unique_id", all.x = TRUE)
  message("Integrated data rows: ", nrow(integrated))
  
  # Check if required columns are present
  required_cols <- c("SC", "cell_type", "secret_ligand")
  for (col in required_cols) {
    if (!col %in% colnames(integrated)) {
      warning("Integrated data missing column: ", col)
    }
  }
  
  return(integrated)
}

# Visualization generation function
# Creates multi-panel PDF with annotation, cell_type, target column, TPS, and filtered plots
generate_visualizations <- function(data, stage, germ_layer, output_prefix, smooth_x_col, smooth_y_col, colorsheet_path, unified_x_range, unified_y_range, target_column, tps_plots = list(), filtered_plots = list()) {
  # Read color sheet
  colorsheet <- read.csv(colorsheet_path)
  my_color <- colorsheet$color
  names(my_color) <- colorsheet$cell_type
  message("Color sheet read successfully")
  
  # Common theme (no legend)
  common_theme <- theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.background = element_rect(fill = "white"),
          legend.position = "none")
  
  # Generate filename with column name, tps_gcv_5ave.dist info
  lambda_info <- ifelse(is.null(tps_plots$lambda), "gcv", paste0("lambda_", tps_plots$lambda, "ave.dist"))
  distance_info <- ifelse(is.null(tps_plots$distance_multiplier), "10ave.dist", paste0(tps_plots$distance_multiplier, "ave.dist"))
  pdf_file <- paste0(output_prefix, germ_layer, "_", target_column, "_2D_visualizations_tps_", lambda_info, "_filtered_", distance_info, ".pdf")
  pdf(pdf_file, width = 15, height = 5)  # Wide format for side-by-side legends
  
  # 1. annotation.2604
  if ("annotation.2604" %in% colnames(data) && "colr.d" %in% colnames(data)) {
    p1 <- ggplot(data, aes_string(x = smooth_x_col, y = smooth_y_col, color = "annotation.2604")) +
      geom_point(size = 3) +
      scale_color_manual(values = setNames(unique(data$colr.d), unique(data$annotation.2604))) +
      coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
      ggtitle(paste(stage, germ_layer, "- annotation.2604")) +
      xlab(smooth_x_col) + ylab(smooth_y_col) + common_theme
    print(p1)
  }
  
  # 2. cell_type
  if ("cell_type" %in% colnames(data)) {
    p2 <- ggplot(data, aes_string(x = smooth_x_col, y = smooth_y_col, color = "cell_type")) +
      geom_point(size = 3) +
      scale_color_manual(values = my_color) +
      coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
      ggtitle(paste(stage, germ_layer, "- cell_type")) +
      xlab(smooth_x_col) + ylab(smooth_y_col) + common_theme
    print(p2)
  }
  
  # 3. Target column (e.g., secret_ligand or other genes)
  if (target_column %in% colnames(data)) {
    p3 <- ggplot(data, aes_string(x = smooth_x_col, y = smooth_y_col, color = target_column)) +
      geom_point(size = 3) +
      scale_color_viridis() +
      coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
      ggtitle(paste(stage, germ_layer, "- ", target_column)) +
      xlab(smooth_x_col) + ylab(smooth_y_col) + common_theme
    print(p3)
  }
  
  # 4. TPS smooth plot
  if (!is.null(tps_plots$plot)) {
    print(tps_plots$plot)
  }
  
  # 5. Filtered smooth plot
  if (!is.null(filtered_plots$plot)) {
    print(filtered_plots$plot)
  }
  
  # Page 6: All legends arranged side by side
  grid::grid.newpage()
  
  # Extract legends
  legends <- list()
  
  if (exists("p1")) {
    p1_with_legend <- p1 + theme(legend.position = "right") + labs(color = "annotation.2604")
    g1 <- ggplotGrob(p1_with_legend)
    legends[[1]] <- g1$grobs[[which(sapply(g1$grobs, function(x) x$name) == "guide-box")]]
  }
  
  if (exists("p2")) {
    p2_with_legend <- p2 + theme(legend.position = "right") + labs(color = "cell_type")
    g2 <- ggplotGrob(p2_with_legend)
    legends[[2]] <- g2$grobs[[which(sapply(g2$grobs, function(x) x$name) == "guide-box")]]
  }
  
  if (exists("p3")) {
    p3_with_legend <- p3 + theme(legend.position = "right") + labs(color = target_column)
    g3 <- ggplotGrob(p3_with_legend)
    legends[[3]] <- g3$grobs[[which(sapply(g3$grobs, function(x) x$name) == "guide-box")]]
  }
  
  # Add TPS and filtered plot legends
  if (!is.null(tps_plots$plot)) {
    tps_with_legend <- tps_plots$plot + theme(legend.position = "right")
    g4 <- ggplotGrob(tps_with_legend)
    legends[[4]] <- g4$grobs[[which(sapply(g4$grobs, function(x) x$name) == "guide-box")]]
  }
  
  if (!is.null(filtered_plots$plot)) {
    filtered_with_legend <- filtered_plots$plot + theme(legend.position = "right")
    g5 <- ggplotGrob(filtered_with_legend)
    legends[[5]] <- g5$grobs[[which(sapply(g5$grobs, function(x) x$name) == "guide-box")]]
  }
  
  # Arrange side by side (horizontal distribution)
  x_positions <- c(0.15, 0.35, 0.55, 0.75, 0.95)  # Multiple positions
  for (i in seq_along(legends)) {
    if (!is.null(legends[[i]])) {
      vp <- grid::viewport(x = x_positions[i], y = 0.5, width = 0.15, height = 0.8)
      grid::pushViewport(vp)
      grid::grid.draw(legends[[i]])
      grid::popViewport()
    }
  }
  
  dev.off()
  message("Visualization PDF generated: ", pdf_file)
  return(NULL)
}


##################################################################




######
######---- Main function: Spatial smoothing for three germ layers ----
main_spatial_smooth <- function(
    input_files = list(),                   # List of 2D coordinate files for three germ layers
    score_file = NULL,                      # Path to 3D_score file
    stage = NULL,                           # Stage identifier
    prefix = NULL,                          # File prefix
    lambda_multiplier = NULL,               # TPS lambda multiplier, default NULL uses GCV auto-selection
    distance_multiplier = NULL,             # Distance filter multiplier, default NULL uses 10x, can be a list
    smooth_x_col = "x_larp",                # Smooth X coordinate column name
    smooth_y_col = "y_larp",                # Smooth Y coordinate column name
    colorsheet_path = "../extend_altas_colorsheet.csv", # Color sheet file path
    target_columns = c("secret_ligand")     # List of target columns to process
) {

  # Parameter validation
  if (length(input_files) == 0) {
    stop("Must provide input_files parameter with 2D coordinate files for three germ layers")
  }
  
  if (is.null(score_file)) {
    stop("Must provide score_file parameter with path to 3D_score file")
  }
  
  if (is.null(stage)) {
    stop("Must provide stage parameter with stage identifier")
  }

  # Generate file prefix
  fn <- prefix
  if (is.null(fn)) {
    fn <- stage
  }
  
  # Create output directory
  dir.create(stage, showWarnings = FALSE)
  setwd(stage)
  message("Output directory: ", getwd())

  if (!file.exists(score_file)) {
    stop("Specified 3D_score file does not exist: ", score_file)
  }
  
  for (file in input_files) {
    if (!file.exists(file)) {
      stop("Specified 2D coordinate file does not exist: ", file)
    }
  }
  # Process each germ layer
  all_data <- list()
  germ_layers <- c("Endoderm", "Mesoderm", "Ectoderm")
  
  # Calculate unified coordinate range across all germ layers first
  message("\n=== Calculating unified coordinate range ===")
  all_x <- c()
  all_y <- c()
  
  for (i in seq_along(input_files)) {
    germ_layer <- germ_layers[i]
    coord_file <- input_files[i]
    
    # Read 2D coordinate file
    coord <- read.csv(coord_file)
    
    # Extract coordinate data
    if (smooth_x_col %in% colnames(coord) && smooth_y_col %in% colnames(coord)) {
      all_x <- c(all_x, coord[[smooth_x_col]])
      all_y <- c(all_y, coord[[smooth_y_col]])
    }
  }
  
  # Calculate unified coordinate range (all germ layers)
  x_min <- min(all_x, na.rm = TRUE)
  x_max <- max(all_x, na.rm = TRUE)
  y_min <- min(all_y, na.rm = TRUE)
  y_max <- max(all_y, na.rm = TRUE)
  
  # Add padding (disabled)
  # x_pad <- (x_max - x_min) * 0.1
  # y_pad <- (y_max - y_min) * 0.1
  unified_x_range <- c(x_min, x_max)
  unified_y_range <- c(y_min, y_max)
  
  message("X: [", round(unified_x_range[1], 2), ", ", round(unified_x_range[2], 2), "]")
  message("Y: [", round(unified_y_range[1], 2), ", ", round(unified_y_range[2], 2), "]")

  # Process each germ layer
  for (i in seq_along(input_files)) {
    germ_layer <- germ_layers[i]
    coord_file <- input_files[i]
    
    message(paste0("\n=== Processing germ layer: ", germ_layer, " ==="))
    
    # Integrate data
    integrated_data <- integrate_data(coord_file, score_file, stage)
    
    # Store data
    all_data[[germ_layer]] <- integrated_data
  }

  # Initialize filtered plot storage
  all_filtered_plots <- list()  # Store filtered plots for all germ layers

  ### ==== Smoothing processing ====
  message("\n=== Starting smoothing processing ===")

  for (i in seq_along(germ_layers)) {
    germ_layer <- germ_layers[i]
    message(paste0("\n=== Processing germ layer: ", germ_layer, " ==="))
    
    # Get data for current germ layer
    base_temp <- all_data[[germ_layer]]
    
    # Validate smooth coordinate columns exist
    if (!(smooth_x_col %in% colnames(base_temp))) {
      stop("Specified smooth X coordinate column '", smooth_x_col, "' does not exist in data")
    }
    if (!(smooth_y_col %in% colnames(base_temp))) {
      stop("Specified smooth Y coordinate column '", smooth_y_col, "' does not exist in data")
    }
    message("Using smooth coordinate columns: X=", smooth_x_col, ", Y=", smooth_y_col)

    # Calculate average nearest neighbor distance (once per germ layer, shared by all columns)
    message("Calculating average nearest neighbor distance...")
    xyrat <- data.frame(x = base_temp[[smooth_x_col]], y = base_temp[[smooth_y_col]])
    nndis <- calc_avg_nnd(xyrat)
    message("Average nearest neighbor distance: ", round(nndis$ave.dist, 3))

    # Determine distance filter multiplier for current germ layer
    if (is.list(distance_multiplier) && length(distance_multiplier) >= i) {
      distance_multiplier_used <- distance_multiplier[[i]]
    } else {
      distance_multiplier_used <- ifelse(is.null(distance_multiplier), 10, distance_multiplier)
    }
    message("Using distance filter multiplier: ", distance_multiplier_used)

    # Create grid for interpolation (once per germ layer, shared by all columns)
    x_range <- range(base_temp[[smooth_x_col]])
    y_range <- range(base_temp[[smooth_y_col]])
    # Use average nearest neighbor distance to determine grid resolution
    grid_resolution <- max(50, round(diff(x_range) / (nndis$ave.dist * 2)))
    x_seq <- seq(x_range[1], x_range[2], length.out = grid_resolution)
    y_seq <- seq(y_range[1], y_range[2], length.out = grid_resolution)
    grid <- expand.grid(x = x_seq, y = y_seq)
    message("Using grid resolution: ", grid_resolution, "x", grid_resolution)

    # Prepare coordinates (once per germ layer, shared by all columns)
    coords <- as.matrix(base_temp[, c(smooth_x_col, smooth_y_col)])

    # Store filtered plots for current germ layer
    germ_layer_filtered_plots <- list()

    # Process each target column
    for (target_col in target_columns) {
      message(paste0("\nProcessing column: ", target_col))
      
      # Validate target column exists
      if (!(target_col %in% colnames(base_temp))) {
        warning("Target column '", target_col, "' does not exist in data, skipping")
        next
      }

      # Prepare values for current column
      values <- base_temp[[target_col]]

      # Initialize variables
      lambda_multiplier_used <- lambda_multiplier
      lambda_used <- NULL
      df_raster <- NULL
      tps_plot <- NULL
      filtered_plot <- NULL

      if (is.null(lambda_multiplier)) {
        # Generate default GCV interpolation only
        message("Performing Thin Plate Spline interpolation (GCV auto-selects lambda)...")
        tps_fit <- Tps(coords, values)
        pred <- predict(tps_fit, grid)

        df_raster <- data.frame(x = grid$x, y = grid$y, density = pred)
        
        # Smooth plot (no legend)
        tps_plot <- ggplot() +
          geom_raster(data = df_raster, aes(x = x, y = y, fill = density), na.rm = FALSE) +
          scale_fill_viridis_c(name = target_col) +
          coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
          ggtitle(paste(stage, germ_layer, "- TPS Interpolation (GCV) -", target_col)) +
          xlab(smooth_x_col) +
          ylab(smooth_y_col) +
          theme_bw() +
          theme(panel.grid = element_blank(), 
                panel.background = element_rect(fill = "white"),
                legend.position = "none")

      } else {
        # Generate manual lambda interpolation
        lambda_used <- nndis$ave.dist * lambda_multiplier_used
        message("Using manual lambda: ", lambda_used, " (multiplier: ", lambda_multiplier_used, ")")

        tps_fit <- Tps(coords, values, lambda = lambda_used)
        pred <- predict(tps_fit, grid)

        df_raster <- data.frame(x = grid$x, y = grid$y, density = pred)
        
        # Smooth plot (no legend)
        tps_plot <- ggplot() +
          geom_raster(data = df_raster, aes(x = x, y = y, fill = density), na.rm = FALSE) +
          scale_fill_viridis_c(name = target_col) +
          coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
          ggtitle(paste(stage, germ_layer, "- TPS Interpolation, lambda =", lambda_used, "-", target_col)) +
          xlab(smooth_x_col) +
          ylab(smooth_y_col) +
          theme_bw() +
          theme(panel.grid = element_blank(), 
                panel.background = element_rect(fill = "white"),
                legend.position = "none")
      }

      ### ==== Part 3: Spatial filtering visualization ====
      message("=== Generating filtered smooth plot ===")

      library(FNN)

      # Calculate nearest neighbor distance for filtering
      original_points <- as.matrix(base_temp[, c(smooth_x_col, smooth_y_col)])
      raster_coords <- as.matrix(df_raster[, c("x", "y")])

      nn_dist <- get.knnx(
        data = original_points,
        query = raster_coords,
        k = 1
      )$nn.dist[, 1]

      # Distance filtering
      distance_threshold <- round(nndis$ave.dist, 3) * distance_multiplier_used
      message("Distance threshold: ", distance_threshold, " (multiplier: ", distance_multiplier_used, ")")

      df_raster_filtered <- df_raster
      df_raster_filtered$density[nn_dist > distance_threshold] <- NA

      # Set title based on whether manual lambda is used
      if (is.null(lambda_multiplier)) {
        title_suffix_filtered <- paste0(target_col, "_near | GCV | threshold=", distance_multiplier_used, "ave.dist")
      } else {
        title_suffix_filtered <- paste0(target_col, "_near | lambda=", lambda_multiplier_used, "ave.dist | threshold=", distance_multiplier_used, "ave.dist")
      }

      # Filtered plot (no legend)
      filtered_plot <- ggplot() +
        geom_raster(data = df_raster_filtered, aes(x = x, y = y, fill = density), na.rm = FALSE) +
        scale_fill_viridis_c(name = target_col, na.value = "transparent") +
        coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
        ggtitle(paste(stage, germ_layer, "- ", title_suffix_filtered)) +
        xlab(smooth_x_col) +
        ylab(smooth_y_col) +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(fill = "white"),
              legend.position = "none")

      # Create small-sized filtered plot for summary
      summary_filtered_plot <- ggplot() +
        geom_raster(data = df_raster_filtered, aes(x = x, y = y, fill = density), na.rm = FALSE) +
        scale_fill_viridis_c(name = target_col, na.value = "transparent") +
        coord_fixed(ratio = 1, xlim = unified_x_range, ylim = unified_y_range) +
        ggtitle(target_col) +
        xlab(smooth_x_col) +
        ylab(smooth_y_col) +
        theme_bw() +
        theme(panel.grid = element_blank(), 
              panel.background = element_rect(fill = "white"),
              legend.position = "none",
              axis.title = element_text(size = 8),
              axis.text = element_text(size = 6),
              plot.title = element_text(size = 10, hjust = 0.5))

      # Store filtered plot
      germ_layer_filtered_plots[[target_col]] <- summary_filtered_plot

      # Generate visualization with all plots
      generate_visualizations(
        data = base_temp,
        stage = stage,
        germ_layer = germ_layer,
        output_prefix = fn,
        smooth_x_col = smooth_x_col,
        smooth_y_col = smooth_y_col,
        colorsheet_path = colorsheet_path,
        unified_x_range = unified_x_range,
        unified_y_range = unified_y_range,
        target_column = target_col,
        tps_plots = list(
          plot = tps_plot,
          lambda = lambda_used,
          distance_multiplier = distance_multiplier_used
        ),
        filtered_plots = list(
          plot = filtered_plot
        )
      )

      message("Column ", target_col, " processing complete!")
    }

    # Save filtered plots for current germ layer
    all_filtered_plots[[germ_layer]] <- germ_layer_filtered_plots

    # Save integrated data
    output_filename <- paste0(fn, "_", germ_layer, "_ST_2D_coord_cbind_mapped_SC_cell_type_meta_score.csv")
    write.csv(base_temp, file = output_filename)
    message("Integrated data saved as: ", output_filename)

    message("Germ layer ", germ_layer, " processing complete!")
  }

  # Generate summary PDF for each germ layer with all target column spatial filtered plots
  for (i in seq_along(germ_layers)) {
    germ_layer <- germ_layers[i]
    message(paste0("\n=== Generating summary PDF for germ layer: ", germ_layer, " ==="))
    
    # Get filtered plots for current germ layer
    filtered_plots_list <- all_filtered_plots[[germ_layer]]
    
    # Determine distance filter multiplier for current germ layer
    if (is.list(distance_multiplier) && length(distance_multiplier) >= i) {
      distance_multiplier_used <- distance_multiplier[[i]]
    } else {
      distance_multiplier_used <- ifelse(is.null(distance_multiplier), 10, distance_multiplier)
    }
    
    # Generate summary PDF
    if (length(filtered_plots_list) > 0) {
      # Calculate layout
      n_plots <- length(filtered_plots_list)
      n_cols <- min(4, ceiling(sqrt(n_plots)))
      n_rows <- ceiling(n_plots / n_cols)
      
      # Generate filename
      lambda_info <- ifelse(is.null(lambda_multiplier), "gcv", paste0("lambda_", lambda_multiplier, "ave.dist"))
      distance_info <- paste0(distance_multiplier_used, "ave.dist")
      pdf_file <- paste0(fn, "_", germ_layer, "_all_score_2D_filtered_summary_tps_", lambda_info, "_filtered_", distance_info, ".pdf")
      
      # Set PDF size
      pdf_width <- max(15, n_cols * 4)
      pdf_height <- max(10, n_rows * 3)
      pdf(pdf_file, width = pdf_width, height = pdf_height)
      
      # Arrange all plots
      gridExtra::grid.arrange(grobs = filtered_plots_list, ncol = n_cols, nrow = n_rows)
      
      dev.off()
      message("Summary PDF generated: ", pdf_file)
    }
  }

  message("\n=== Analysis complete! ===")
  message("Generated files:")
  message(" - Integrated data files for each germ layer")
  message(" - Visualization PDF files for each target column in each germ layer")
  message(" - Summary PDF files with all target column filtered plots for each germ layer")

  # Switch back to original working directory
  setwd("..")
  message("Switched back to original working directory: ", getwd())

  return(list(
    all_data = all_data,
    unified_x_range = unified_x_range,
    unified_y_range = unified_y_range
  ))
}

# ===== Usage Examples =====

# ⭐ Example 1 (Recommended): Process three germ layers with same distance filter multiplier
# result1 <- main_spatial_smooth(
#   input_files = list(
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Endoderm_LARP.2d.coord.csv",
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Mesoderm_Midline_LARP.2d.coord.csv",
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Ectoderm_LARP.2d.coord.csv"
#   ),
#   score_file = "E7.25_3D_score.csv",
#   stage = "E7.25",
#   prefix = "E7.25_mapped_sc_secrete_gene_set_AUC",
#   # lambda_multiplier = NULL,             # Default NULL, uses GCV auto-selection
#   distance_multiplier = 5
# )

# ⭐ Example 2: Process three germ layers with different distance filter multipliers
# result3 <- main_spatial_smooth(
#   input_files = list(
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Endoderm_LARP.2d.coord.csv",
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Mesoderm_Midline_LARP.2d.coord.csv",
#     "E7.25_mapped_sc_secrete_gene_set_AUC_Ectoderm_LARP.2d.coord.csv"
#   ),
#   score_file = "E7.25_3D_score.csv",
#   stage = "E7.25",
#   prefix = "E7.25_mapped_sc_secrete_gene_set_AUC",
#   # lambda_multiplier = NULL,             # Default NULL, uses GCV auto-selection
#   distance_multiplier = list(5, 8, 10)  # Specify different distance filter multipliers for each germ layer
# )

# ⭐ Example 3: Use different coordinate columns
# result2 <- main_spatial_smooth(
#   input_files = list(
#     "E7.5_mapped_sc_secrete_gene_set_AUC_Endoderm_with_2d_coord.csv",
#     "E7.5_mapped_sc_secrete_gene_set_AUC_Mesoderm_Midline_with_2d_coord.csv",
#     "E7.5_mapped_sc_secrete_gene_set_AUC_Ectoderm_with_2d_coord.csv"
#   ),
#   score_file = "E7.5_3D_score.csv",
#   stage = "E7.5",
#   prefix = "E7.5_mapped_sc_secrete_gene_set_AUC",
#   distance_multiplier = list(6, 9, 12),  # Specify different distance filter multipliers for each germ layer
#   smooth_x_col="x2d",
#   smooth_y_col="y2d"
# )


