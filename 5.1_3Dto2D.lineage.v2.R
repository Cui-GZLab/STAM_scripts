#!/usr/bin/env Rscript

## 3Dto2D.lineage.v2.R

# ==============================================================================
# 3D to 2D Projection for Lineage Analysis
# ==============================================================================
#
# Description:
# This script projects 3D spatial coordinates onto a 2D plane based on the
# Anterior-Posterior (A-P) axis. It calculates the A-P axis slope from centroid
# positions, projects each sample point onto the axis line, and generates
# rotated coordinates where the A-P axis is horizontal. 
#
# Input files:
#   - RDS file: 3D coordinates (x, y, z) for each sample (rand.Real_points2mesh.rds)
#   - CSV file: Metadata containing stage, sectorlcm, lineagelcm, condition,
#               colr.s (sector colors), colr.d (domain colors), lineage.plot
#
# Output files:
#   - CSV files: 2D projected coordinates for each lineage (*_xy.csv)
#   - PNG files: Visualization plots showing projected points colored by
#               sectorlcm or annotation (*_2d.sectorlcm.line.png, *_larp.*.png)
#   - HTML file: 3D visualization of A-P axis sectors (AP_sector.html)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(threejs)
  library(htmlwidgets)
  library(dplyr)
  library(plotly)
  library(ggplot2)
})

# ============================================================================
# Helper Functions
# ============================================================================

#' Calculate two points on a line at distance z from a given point
#' 
#' @param x0 Known point x coordinate
#' @param y0 Known point y coordinate
#' @param z Distance from the point
#' @param k Line slope
#' @return Returns four values: p1_x, p1_y, p2_x, p2_y
get_two_points_on_line <- function(x0, y0, z, k) {
  if (!is.numeric(x0) || !is.numeric(y0)) {
    stop("Error: Known point (x0, y0) must be numeric")
  }
  if (!is.numeric(k) && !is.infinite(k)) {
    stop("Error: Slope k must be numeric or Inf")
  }
  if (!is.numeric(z) || z <= 0) {
    stop("Error: Distance z must be positive")
  }
  
  if (is.infinite(k)) {
    p1_x <- x0
    p1_y <- y0 + z
    p2_x <- x0
    p2_y <- y0 - z
  } else {
    norm <- sqrt(1 + k^2)
    ux <- 1 / norm
    uy <- k / norm
    p1_x <- x0 + z * ux
    p1_y <- y0 + z * uy
    p2_x <- x0 - z * ux
    p2_y <- y0 - z * uy
  }
  return(c(p1_x = p1_x, p1_y = p1_y, p2_x = p2_x, p2_y = p2_y))
}

#' Determine point position relative to a line segment (left or right)
#' 
#' @param point_x Point x coordinate
#' @param point_y Point y coordinate
#' @param line_start_x Line segment start x coordinate
#' @param line_start_y Line segment start y coordinate
#' @param line_end_x Line segment end x coordinate
#' @param line_end_y Line segment end y coordinate
#' @return 1 for left side, -1 for right side, 0 for on the line
#'


#' Determine the geometric left/right position of a point relative to AP axis
#' 
#' @return 1 for left side of AP axis, -1 for right side, 0 for on the axis
get_point_side_geometric <- function(point_x, point_y, 
                                     line_start_x, line_start_y, 
                                     line_end_x, line_end_y) {
  
  # Calculate midpoint of AP axis
  mid_x <- (line_start_x + line_end_x) / 2
  mid_y <- (line_start_y + line_end_y) / 2
  
  # Calculate slope of AP axis
  ap_slope <- (line_end_y - line_start_y) / (line_end_x - line_start_x)
  
  # Calculate perpendicular bisector slope (perpendicular to AP axis)
  if (ap_slope == 0) {
    # AP axis is horizontal, perpendicular bisector is vertical
    # Perpendicular bisector equation: x = mid_x
    # Compare point's x coordinate with mid_x
    if (point_x < mid_x) {
      return(1)  # Left side of perpendicular bisector
    } else if (point_x > mid_x) {
      return(-1) # Right side of perpendicular bisector
    } else {
      return(0)  # On the perpendicular bisector
    }
  } else if (is.infinite(ap_slope)) {
    # AP axis is vertical, perpendicular bisector is horizontal
    # Perpendicular bisector equation: y = mid_y
    # Compare point's y coordinate with mid_y
    # Note: Need to determine left/right definition
    # Assume AP axis is vertical, upward is positive direction
    # Left side: x < mid_x
    if (point_x < mid_x) {
      return(1)  # Left side of perpendicular bisector
    } else if (point_x > mid_x) {
      return(-1) # Right side of perpendicular bisector
    } else {
      return(0)  # On the perpendicular bisector
    }
  } else {
    # Perpendicular bisector slope
    perpendicular_slope <- -1 / ap_slope
    
    # Perpendicular bisector equation: y - mid_y = perpendicular_slope * (x - mid_x)
    # For target point (point_x, point_y), find point on perpendicular bisector with same y
    # Solve: point_y - mid_y = perpendicular_slope * (x - mid_x)
    # Get: x = mid_x + (point_y - mid_y) / perpendicular_slope
    x_on_perpendicular <- mid_x + (point_y - mid_y) / perpendicular_slope
    
    # Compare target point's x coordinate with perpendicular bisector point
    if (point_x < x_on_perpendicular) {
      return(1)  # Left side of perpendicular bisector
    } else if (point_x > x_on_perpendicular) {
      return(-1) # Right side of perpendicular bisector
    } else {
      return(0)  # On the perpendicular bisector
    }
  }
}
#' Calculate A-P axis slope
#' 
#' @param sap3d 3D data frame
#' @param ap_sectors A-P axis sector names, e.g., c("A", "P")
#' @param sector_col Sector column name
#' @param condition_col Condition column name
#' @return Returns slope K and centroid coordinates
#'
calculate_ap_slope <- function(sap3d, ap_sectors = c("A", "P"), 
                                sector_col = "sectorlcm", 
                                condition_col = "condition") {
  # Filter A-P region
  sap3df <- filter(sap3d, .data[[sector_col]] %in% ap_sectors)
  
  if (nrow(sap3df) == 0) {
    stop(paste("No A-P region data found, please check if", sector_col, "column contains", paste(ap_sectors, collapse = ", ")))
  }
  
  # Step 1: Calculate centroids by condition
  centroid_list <- list()
  for (s in unique(sap3df[[condition_col]])) {
    fil1 <- filter(sap3df, .data[[condition_col]] == s)
    centroid_xy <- apply(fil1[, c("x", "y", "z")], 2, mean)
    name_val <- unique(fil1$name)[1]
    centroid_list[[name_val]] <- centroid_xy
  }
  
  centroid_df <- do.call(rbind, centroid_list)
  colnames(centroid_df) <- c("x", "y", "z")
  
  # Add sectorlcm information
  pho1 <- unique(sap3df[, c("name", sector_col)])
  centroid_df <- merge(
    as.data.frame(centroid_df),
    pho1,
    by.x = 0,
    by.y = "name",
    all.x = TRUE
  )
  
  # Step 2: Calculate final centroids for A and P by averaging again
  xy0 <- data.frame(x = NULL, y = NULL)
  for (s in unique(centroid_df[[sector_col]])) {
    fil1 <- filter(centroid_df, .data[[sector_col]] == s)
    centroid_xy <- apply(fil1[, c("x", "y")], 2, mean)
    centroid_xy_t <- as.data.frame(t(centroid_xy))
    rownames(centroid_xy_t) <- s
    xy0 <- rbind(xy0, centroid_xy_t)
  }
  colnames(xy0) <- c("x", "y")
  
  # Ensure both A and P exist
  if (!all(ap_sectors %in% rownames(xy0))) {
    stop(paste("Missing A or P region centroid data"))
  }
  
  # Calculate slope (using xy0)
  K <- (xy0[ap_sectors[1], "y"] - xy0[ap_sectors[2], "y"]) / 
       (xy0[ap_sectors[1], "x"] - xy0[ap_sectors[2], "x"])
  
  message("A-P axis slope K = ", round(K, 4))
  
  return(list(
    K = K,
    centroids = xy0,
    ap_data = sap3df
  ))
}

#' Project 3D coordinates to 2D plane
#' 
#' @param sap3d 3D data frame
#' @param lineage_filter Lineage value to filter
#' @param K A-P axis slope
#' @param centroids A-P axis centroid coordinates (contains A and P coordinates)
#' @param anterior_sector Anterior sector identifier
#' @param posterior_sector Posterior sector identifier
#' @param use_perpendicular Whether to use perpendicular bisector for direction (TRUE: use perpendicular, FALSE: use sector name)
#' @param ap_sectors A-P axis sector identifiers, default c("A", "P")
#' @param sector_col Sector column name
#' @param lineage_col Lineage column name
#' @return Returns data frame with 2D coordinates
#'
project_3d_to_2d <- function(sap3d, lineage_filter, K, 
                              centroids = NULL,
                              anterior_sector = NULL,
                              posterior_sector = NULL,
                              use_perpendicular = TRUE,
                              ap_sectors = c("A", "P"),
                              sector_col = "sectorlcm",
                              lineage_col = "lineagelcm") {
  
  # Use passed sap3d directly, no germ layer filtering
  # Because filtering was already done based on lineage.plot column in main function
  sap3df <- sap3d
  
  if (nrow(sap3df) == 0) {
    warning("No data found")
    return(NULL)
  }
  
  message("Processing data points: ", nrow(sap3df))
  
  # Prepare projection data
  dfxyz <- sap3df[, c("x", "y", "z")]
  dfxyz$k <- K
  
  # Execute projection
  result_matrix <- apply(
    X = dfxyz,
    MARGIN = 1,
    FUN = function(row) {
      get_two_points_on_line(
        x0 = row[1],
        y0 = row[2],
        z = row[3],
        k = row[4]
      )
    }
  )
  
  xyap <- as.data.frame(t(result_matrix))
  colnames(xyap) <- c("p1_x", "p1_y", "p2_x", "p2_y")
  sap3df <- cbind(sap3df, xyap)
  
  # Show sector classification information
  sectors_unique <- unique(sap3df[[sector_col]])
  message("    Detected sectors: ", paste(sectors_unique, collapse = ", "))
  
  # Determine direction based on parameters
  if (!is.null(anterior_sector) && !is.null(posterior_sector)) {
    # Method 1: Use specified anterior/posterior sectors
    message("    Using specified sectors to determine direction")
    sap3df$x2d <- ifelse(sap3df[[sector_col]] == anterior_sector, sap3df$p2_x, sap3df$p1_x)
    sap3df$y2d <- ifelse(sap3df[[sector_col]] == anterior_sector, sap3df$p2_y, sap3df$p1_y)
  } else if (use_perpendicular && !is.null(centroids) && nrow(centroids) >= 2) {
    # Method 2: Use perpendicular bisector
    # Get A and P coordinates
    ap_points <- centroids[1:2, ]
    line_start_x <- ap_points[1, "x"]
    line_start_y <- ap_points[1, "y"]
    line_end_x <- ap_points[2, "x"]
    line_end_y <- ap_points[2, "y"]
    
    message("    Using perpendicular bisector to determine direction")
    message("    A-P axis start: (", round(line_start_x, 2), ", ", round(line_start_y, 2), ")")
    message("    A-P axis end: (", round(line_end_x, 2), ", ", round(line_end_y, 2), ")")
    
    # Calculate geometric left/right position for each point relative to AP axis
    sap3df$side <- mapply(
      get_point_side_geometric,
      point_x = sap3df$x,
      point_y = sap3df$y,
      line_start_x = line_start_x,
      line_start_y = line_start_y,
      line_end_x = line_end_x,
      line_end_y = line_end_y
    )
    
    # Count points on each side
    left_count <- sum(sap3df$side == 1)
    right_count <- sum(sap3df$side == -1)
    on_line_count <- sum(sap3df$side == 0)
    message("    Points on one side of perpendicular bisector: ", left_count)
    message("    Points on other side of perpendicular bisector: ", right_count)
    message("    Points on perpendicular bisector: ", on_line_count)
    
    # Select projection point based on position: left side uses p2, right side uses p1
    sap3df$x2d <- ifelse(sap3df$side == 1, sap3df$p2_x, sap3df$p1_x)
    sap3df$y2d <- ifelse(sap3df$side == 1, sap3df$p2_y, sap3df$p1_y)
    
    # Rotation: Rotate AP axis to horizontal, A on left, P on right
    if (!is.null(centroids) && nrow(centroids) >= 2) {
      # Get A and P coordinates
      a_row <- which(rownames(centroids) %in% ap_sectors[1])
      p_row <- which(rownames(centroids) %in% ap_sectors[2])
      
      if (length(a_row) > 0 && length(p_row) > 0) {
        a_x <- centroids[a_row, "x"]
        a_y <- centroids[a_row, "y"]
        p_x <- centroids[p_row, "x"]
        p_y <- centroids[p_row, "y"]
        
        # Calculate AP vector angle (relative to x axis)
        dx <- p_x - a_x
        dy <- p_y - a_y
        angle <- atan2(dy, dx)
        
        # Rotation matrix: Rotate clockwise by angle to make AP axis horizontal
        cos_theta <- cos(-angle)
        sin_theta <- sin(-angle)
        
        # Rotate all points, using A as origin
        sap3df$x_larp <- (sap3df$x2d - a_x) * cos_theta - (sap3df$y2d - a_y) * sin_theta
        sap3df$y_larp <- (sap3df$x2d - a_x) * sin_theta + (sap3df$y2d - a_y) * cos_theta
        
        message("    Rotation angle: ", round(angle * 180 / pi, 2), " degrees")
        message("    Rotated A point coordinates: (0, 0)")
        message("    Rotated P point coordinates: (", round(sqrt(dx^2 + dy^2), 2), ", 0)")
      }
    }
  } else {
    # Method 3: Determine based on sector name containing A
    message("    Using sector name to determine direction")
    
    # Determine left/right position of A and P centroids
    if (!is.null(centroids) && nrow(centroids) >= 2) {
      # Get A and P coordinates
      a_row <- which(rownames(centroids) %in% ap_sectors[1])
      p_row <- which(rownames(centroids) %in% ap_sectors[2])
      
      if (length(a_row) > 0 && length(p_row) > 0) {
        a_x <- centroids[a_row, "x"]
        p_x <- centroids[p_row, "x"]
        
        if (a_x < p_x) {
          message("    A centroid on left (x = ", round(a_x, 2), "), P centroid on right (x = ", round(p_x, 2), ")")
          # A on left, sectors containing A use p2 direction
          sectors_with_A <- sectors_unique[grepl("A", sectors_unique)]
          sectors_without_A <- sectors_unique[!grepl("A", sectors_unique)]
          message("    Sectors with A (use p2): ", paste(sectors_with_A, collapse = ", "))
          message("    Sectors without A (use p1): ", paste(sectors_without_A, collapse = ", "))
          sap3df$x2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p2_x, sap3df$p1_x)
          sap3df$y2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p2_y, sap3df$p1_y)
        } else {
          message("    P centroid on left (x = ", round(p_x, 2), "), A centroid on right (x = ", round(a_x, 2), ")")
          # P on left, sectors containing A use p1 direction
          sectors_with_A <- sectors_unique[grepl("A", sectors_unique)]
          sectors_without_A <- sectors_unique[!grepl("A", sectors_unique)]
          message("    Sectors with A (use p1): ", paste(sectors_with_A, collapse = ", "))
          message("    Sectors without A (use p2): ", paste(sectors_without_A, collapse = ", "))
          sap3df$x2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_x, sap3df$p2_x)
          sap3df$y2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_y, sap3df$p2_y)
        }
      } else {
        # Default logic: sectors with A use p1, without A use p2
        sectors_with_A <- sectors_unique[grepl("A", sectors_unique)]
        sectors_without_A <- sectors_unique[!grepl("A", sectors_unique)]
        message("    Sectors with A (use p1): ", paste(sectors_with_A, collapse = ", "))
        message("    Sectors without A (use p2): ", paste(sectors_without_A, collapse = ", "))
        sap3df$x2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_x, sap3df$p2_x)
        sap3df$y2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_y, sap3df$p2_y)
      }
      
      # Rotation: Rotate AP axis to horizontal, A on left, P on right
      if (!is.null(centroids) && nrow(centroids) >= 2) {
        # Get A and P coordinates
        a_row <- which(rownames(centroids) %in% ap_sectors[1])
        p_row <- which(rownames(centroids) %in% ap_sectors[2])
        
        if (length(a_row) > 0 && length(p_row) > 0) {
          a_x <- centroids[a_row, "x"]
          a_y <- centroids[a_row, "y"]
          p_x <- centroids[p_row, "x"]
          p_y <- centroids[p_row, "y"]
          
          # Calculate AP vector angle (relative to x axis)
          dx <- p_x - a_x
          dy <- p_y - a_y
          angle <- atan2(dy, dx)
          
          # Rotation matrix: Rotate clockwise by angle to make AP axis horizontal
          cos_theta <- cos(-angle)
          sin_theta <- sin(-angle)
          
          # Rotate all points, using A as origin
          sap3df$x_larp <- (sap3df$x2d - a_x) * cos_theta - (sap3df$y2d - a_y) * sin_theta
          sap3df$y_larp <- (sap3df$x2d - a_x) * sin_theta + (sap3df$y2d - a_y) * cos_theta
          
          message("    Rotation angle: ", round(angle * 180 / pi, 2), " degrees")
          message("    Rotated A point coordinates: (0, 0)")
          message("    Rotated P point coordinates: (", round(sqrt(dx^2 + dy^2), 2), ", 0)")
        }
      }
    } else {
      # Default logic: sectors with A use p1, without A use p2
      sectors_with_A <- sectors_unique[grepl("A", sectors_unique)]
      sectors_without_A <- sectors_unique[!grepl("A", sectors_unique)]
      message("    Sectors with A (use p1): ", paste(sectors_with_A, collapse = ", "))
      message("    Sectors without A (use p2): ", paste(sectors_without_A, collapse = ", "))
      sap3df$x2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_x, sap3df$p2_x)
      sap3df$y2d <- ifelse(grepl("A", sap3df[[sector_col]]), sap3df$p1_y, sap3df$p2_y)
      
      # Rotation: Rotate AP axis to horizontal, A on left, P on right
      if (!is.null(centroids) && nrow(centroids) >= 2) {
        # Get A and P coordinates
        a_row <- which(rownames(centroids) %in% ap_sectors[1])
        p_row <- which(rownames(centroids) %in% ap_sectors[2])
        
        if (length(a_row) > 0 && length(p_row) > 0) {
          a_x <- centroids[a_row, "x"]
          a_y <- centroids[a_row, "y"]
          p_x <- centroids[p_row, "x"]
          p_y <- centroids[p_row, "y"]
          
          # Calculate AP vector angle (relative to x axis)
          dx <- p_x - a_x
          dy <- p_y - a_y
          angle <- atan2(dy, dx)
          
          # Rotation matrix: Rotate clockwise by angle to make AP axis horizontal
          cos_theta <- cos(-angle)
          sin_theta <- sin(-angle)
          
          # Rotate all points, using A as origin
          sap3df$x_larp <- (sap3df$x2d - a_x) * cos_theta - (sap3df$y2d - a_y) * sin_theta
          sap3df$y_larp <- (sap3df$x2d - a_x) * sin_theta + (sap3df$y2d - a_y) * cos_theta
          
          message("    Rotation angle: ", round(angle * 180 / pi, 2), " degrees")
          message("    Rotated A point coordinates: (0, 0)")
          message("    Rotated P point coordinates: (", round(sqrt(dx^2 + dy^2), 2), ", 0)")
        }
      }
    }
  }
  
  return(sap3df)
}

#' Save result files
#' 
#' @param data Data frame
#' @param fn File name prefix
#' @param lineage_name Lineage name (for file name)
#' @param save_html Whether to save visualization
#' @param html_template HTML template data (centroid coordinates)
#' @param x_range X axis range
#' @param y_range Y axis range
#' @param centroids A-P axis centroid coordinates
#' @param ap_sectors A-P axis sector identifiers
#'
save_results <- function(data, fn, lineage_name, save_html = TRUE, html_template = NULL, x_range = NULL, y_range = NULL, centroids = NULL, ap_sectors = c("A", "P")) {
  # Save CSV
  csv_file <- paste0(fn, lineage_name, "_xy.csv")
  write.csv(data, file = csv_file, row.names = FALSE)
  message("  Saved CSV: ", csv_file)
  

  
  # Save visualization (using ggplot2)
  if (save_html && !is.null(html_template)) {
    # Ensure html_template has correct column names
    if (!all(c("x", "y") %in% colnames(html_template))) {
      warning("html_template missing x or y columns, skipping visualization")
    } else if (nrow(html_template) < 2) {
      warning("html_template has insufficient data, skipping visualization")
    } else if (nrow(data) == 0) {
      warning("Data is empty, skipping visualization")
    } else {
      # Ensure ggplot2 package is installed
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        install.packages("ggplot2")
      }
      library(ggplot2)
      
      # 1. Visualize x2d, y2d (colored by sectorlcm, using colr.s color codes)
      if ("colr.s" %in% colnames(data) && "sectorlcm" %in% colnames(data)) {
        p1 <- ggplot() +
          # Add AP axis line
          geom_line(data = html_template, aes(x = x, y = y), color = "#1f77b4", linewidth = 1.5) +
          # Add data points (using colr.s as color)
          geom_point(data = data, aes(x = x2d, y = y2d, color = sectorlcm), size = 3) +
          # Set colors
          scale_color_manual(values = setNames(unique(data$colr.s), unique(data$sectorlcm))) +
          # Set axis ratio to 1:1
          coord_fixed(ratio = 1) +
          # Add title and labels
          ggtitle(paste(lineage_name, " - 2D projection (x2d, y2d) - sectorlcm")) +
          xlab("x") +
          ylab("y") +
          theme_academic()
        
        # Save as PNG
        png_file1 <- paste0(fn, lineage_name, "_2d.sectorlcm.line.png")
        ggsave(png_file1, plot = p1, width = 10, height = 2.5, dpi = 300)
        message("  Saved PNG: ", png_file1)
      }
      
      # 2. Visualize x2d, y2d (colored by annotation.2604, using colr.d color codes)
      if ("colr.d" %in% colnames(data) && "annotation.2604" %in% colnames(data)) {
        p2 <- ggplot() +
          # Add AP axis line
          geom_line(data = html_template, aes(x = x, y = y), color = "#1f77b4", linewidth = 1.5) +
          # Add data points (using colr.d as color)
          geom_point(data = data, aes(x = x2d, y = y2d, color = annotation.2604), size = 3) +
          # Set colors
          scale_color_manual(values = setNames(unique(data$colr.d), unique(data$annotation.2604))) +
          # Set axis ratio to 1:1
          coord_fixed(ratio = 1) +
          # Add title and labels
          ggtitle(paste(lineage_name, " - 2D projection (x2d, y2d) - annotation.2604")) +
          xlab("x") +
          ylab("y") +
          theme_academic()
        
        # Save as PNG
        png_file2 <- paste0(fn, lineage_name, "_2d.annotation.line.png")
        ggsave(png_file2, plot = p2, width = 10, height = 2.5, dpi = 300)
        message("  Saved PNG: ", png_file2)
      }
      
      # 3. Visualize x_larp, y_larp (if exists)
      if ("x_larp" %in% colnames(data) && "y_larp" %in% colnames(data)) {
        # Calculate rotated AP axis (horizontal)
        if (!is.null(centroids) && nrow(centroids) >= 2) {
          # Get A and P coordinates
          a_row <- which(rownames(centroids) %in% ap_sectors[1])
          p_row <- which(rownames(centroids) %in% ap_sectors[2])
          
          if (length(a_row) > 0 && length(p_row) > 0) {
            a_x <- centroids[a_row, "x"]
            a_y <- centroids[a_row, "y"]
            p_x <- centroids[p_row, "x"]
            p_y <- centroids[p_row, "y"]
            
            # Calculate AP vector angle
            dx <- p_x - a_x
            dy <- p_y - a_y
            angle <- atan2(dy, dx)
            
            # Rotation matrix
            cos_theta <- cos(-angle)
            sin_theta <- sin(-angle)
            
            # Rotate AP axis points
            rotated_ap <- data.frame(
              x = c(0, sqrt(dx^2 + dy^2)),
              y = c(0, 0)
            )
            
            # 1. Colored by sectorlcm (using colr.s color codes)
            if ("colr.s" %in% colnames(data) && "sectorlcm" %in% colnames(data)) {
              p3 <- ggplot() +
                # Add rotated AP axis line
                geom_line(data = rotated_ap, aes(x = x, y = y), color = "#1f77b4", linewidth = 1.5) +
                # Add data points (using colr.s as color)
                geom_point(data = data, aes(x = x_larp, y = y_larp, color = sectorlcm), size = 3) +
                # Set colors
                scale_color_manual(values = setNames(unique(data$colr.s), unique(data$sectorlcm))) +
                # Set axis ratio to 1:1
                coord_fixed(ratio = 1) +
                # Add title and labels
                ggtitle(paste(lineage_name, " - Rotated (x_larp, y_larp) - sectorlcm")) +
                xlab("x (rotated)") +
                ylab("y (rotated)") +
                theme_academic()
              
              # Save as PNG
              png_file3 <- paste0(fn, lineage_name, "_larp.sectorlcm.line.png")
              ggsave(png_file3, plot = p3, width = 10, height = 2.5, dpi = 300)
              message("  Saved PNG: ", png_file3)
            }
            
            # 2. Colored by annotation.2604 (using colr.d color codes)
            if ("colr.d" %in% colnames(data) && "annotation.2604" %in% colnames(data)) {
              p4 <- ggplot() +
                # Add rotated AP axis line
                geom_line(data = rotated_ap, aes(x = x, y = y), color = "#1f77b4", linewidth = 1.5) +
                # Add data points (using colr.d as color)
                geom_point(data = data, aes(x = x_larp, y = y_larp, color = annotation.2604), size = 3) +
                # Set colors
                scale_color_manual(values = setNames(unique(data$colr.d), unique(data$annotation.2604))) +
                # Set axis ratio to 1:1
                coord_fixed(ratio = 1) +
                # Add title and labels
                ggtitle(paste(lineage_name, " - Rotated (x_larp, y_larp) - annotation.2604")) +
                xlab("x (rotated)") +
                ylab("y (rotated)") +
                theme_academic()
              
              # Save as PNG
              png_file4 <- paste0(fn, lineage_name, "_larp.annotation.line.png")
              ggsave(png_file4, plot = p4, width = 10, height = 2.5, dpi = 300)
              message("  Saved PNG: ", png_file4)
            }
          }
        }
      }
    }
  }
}

# ============================================================================
# Main Function
# ============================================================================

#' 3D to 2D projection main function
#' 
#' @param input_rds Input RDS file path (3D coordinates), if NULL will auto-find based on stage
#' @param input_meta Input metadata CSV file path
#' @param stage Stage identifier (e.g., E6.75, E7.0, E7.25, E7.5), used to filter from merged meta file
#' @param rds_dir RDS file directory (used when input_rds is NULL)
#' @param rds_pattern RDS file name matching pattern (default "rand.Real_points2mesh.rds")
#' @param output_dir Output directory
#' @param file_prefix File prefix
#' @param ap_sectors A-P axis sector identifiers, default c("A", "P")
#' @param lineages List of lineages to process, default list(Endoderm = "Endoderm", Mesoderm_Midline = c("Mesoderm", "Midline"))
#' @param use_perpendicular Whether to use perpendicular bisector for direction (TRUE: use perpendicular, FALSE: use sector name)
#' @param sector_col Sector column name
#' @param lineage_col Lineage column name
#' @param condition_col Condition column name
#' @param save_html Whether to save HTML visualization
#' @param lineage_methods Direction determination method for each lineage, format list("Endoderm" = "perpendicular", "Mesoderm" = "sector")
#' @return Returns processing result list
#'
convert_3dto2d <- function(
    input_rds = NULL,
    input_meta,
    stage = NULL,
    rds_dir = NULL,
    rds_pattern = "rand.Real_points2mesh.rds",
    output_dir = NULL,
    file_prefix = NULL,
    ap_sectors = c("A", "P"),
    lineages = list(
      Endoderm = "Endoderm",
      Mesoderm_Midline = c("Mesoderm", "Midline")
    ),
    use_perpendicular = TRUE,
    lineage_methods = NULL,
    sector_col = "sectorlcm",
    lineage_col = "lineagelcm",
    condition_col = "condition",
    save_html = TRUE
) {
  
  # Stage mapping (E7.5 special handling)
  stage_map <- c(
    "E6.75" = "E6.75",
    "E7.0" = "E7.0",
    "E7.25" = "E7.25",
    "E7.5" = "E7.53"
  )
  
  # Check required parameters
  if (!file.exists(input_meta)) {
    stop("Input metadata file does not exist: ", input_meta)
  }
  
  # Set stage and file_prefix
  if (is.null(stage)) {
    if (is.null(file_prefix)) {
      stop("Must provide stage or file_prefix parameter")
    }
    stage <- file_prefix
  }
  if (is.null(file_prefix)) {
    file_prefix <- stage
  }
  
  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- file_prefix
  }
  
  fn0 <- paste0(file_prefix, "_3Dto2d.")
  
  # Save original working directory
  original_wd <- getwd()
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # Get full path of output directory
  output_dir_full <- normalizePath(output_dir)
  setwd(output_dir)
  message("Output directory: ", output_dir_full)
  
  # Read data
  message("\n=== Reading data ===")
  
  # If input_rds not provided, auto-find based on stage (before switching directory)
  if (is.null(input_rds)) {
    if (is.null(rds_dir)) {
      # Use mapped directory name
      actual_stage <- stage_map[stage]
      if (is.na(actual_stage)) {
        actual_stage <- stage
      }
      # Build path based on original working directory
      rds_dir <- file.path(original_wd, "..", actual_stage)
    }
    
    fls <- list.files(
      path = rds_dir,
      pattern = rds_pattern,
      recursive = FALSE,
      full.names = TRUE
    )
    
    if (length(fls) == 0) {
      stop(paste("No RDS file matching pattern", rds_pattern, "found in directory", rds_dir))
    }
    
    input_rds <- fls[1]
    message("Auto-found RDS file: ", input_rds)
  }
  
  if (!file.exists(input_rds)) {
    stop("Input RDS file does not exist: ", input_rds)
  }
  
  message("Reading 3D coordinates: ", input_rds)
  sap3d <- readRDS(input_rds)
  message("3D data points: ", nrow(sap3d))
  
  # Add stage and stname columns
  sap3d$stage <- stage
  sap3d$stname <- paste(sap3d$stage, sap3d$name, sep = "_")
  
  # Read meta file based on original working directory
  meta_path <- file.path(original_wd, input_meta)
  message("Reading metadata: ", meta_path)
  metanew <- read.csv(meta_path)
  message("Metadata total rows: ", nrow(metanew))
  
  # Filter metadata by stage
  if ("stage" %in% colnames(metanew)) {
    metanew <- metanew %>% filter(stage == !!stage)
    message("Filtered metadata rows (stage = ", stage, "): ", nrow(metanew))
  } else {
    warning("No stage column in metadata, using all data")
  }
  
  # Merge data (using stname)
  sap3d <- merge(sap3d, metanew, by = "stname")
  message("Merged data points: ", nrow(sap3d))
  
  # Process lineage.plot column
  if ("lineage.plot" %in% colnames(sap3d)) {
    # Filter out samples with lineage.plot = "remove"
    sap3d <- sap3d %>% filter(lineage.plot != "remove")
    message("Filtered data points: ", nrow(sap3d))
  } else {
    warning("No lineage.plot column in metadata, skipping filter")
  }
  
  # Check required columns
  required_cols <- c("x", "y", "z", "name", sector_col, lineage_col, condition_col)
  missing_cols <- setdiff(required_cols, colnames(sap3d))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Calculate A-P axis
  message("\n=== Calculating A-P axis ===")
  ap_result <- calculate_ap_slope(sap3d, ap_sectors, sector_col, condition_col)
  K <- ap_result$K
  
  # Save A-P axis visualization
  if (save_html) {
    fn_ap <- paste0(fn0, "AP_")
    colrs <- unique(ap_result$ap_data$colr.s)
    names(colrs) <- unique(ap_result$ap_data[[sector_col]])
    
    fig <- plot_ly(ap_result$ap_data, x = ~x, y = ~y, z = ~z, 
                   color = ~sectorlcm, colors = colrs, size = 2)
    saveWidget(fig, file = paste0(fn_ap, "sector.html"), selfcontained = FALSE)
    message("Saved A-P axis visualization: ", paste0(fn_ap, "sector.html"))
  }
  
  # Process each lineage
  message("\n=== Processing lineages ===")
  results <- list()
  all_data <- data.frame()
  
  # Check if lineage.plot column exists
  if ("lineage.plot" %in% colnames(sap3d)) {
    # Get all unique lineage names from lineage.plot column (handling semicolon-separated cases)
    all_lineage_names <- unique(unlist(strsplit(sap3d$lineage.plot, ";")))
    all_lineage_names <- all_lineage_names[all_lineage_names != "remove"]
    
    # Step 1: Process all lineages and collect data
    for (lineage_name in all_lineage_names) {
      # Filter samples belonging to current lineage
      lineage_filter <- lineage_name
      sap3d_subset <- sap3d[grepl(lineage_name, sap3d$lineage.plot, fixed = TRUE), ]
      
      if (nrow(sap3d_subset) == 0) {
        message("\nSkipping lineage: ", lineage_name, " (no data)")
        next
      }
      
      # Determine method for current lineage
      current_method <- use_perpendicular
      if (!is.null(lineage_methods) && lineage_name %in% names(lineage_methods)) {
        method_val <- lineage_methods[[lineage_name]]
        if (method_val == "sector") {
          current_method <- FALSE
        } else if (method_val == "perpendicular") {
          current_method <- TRUE
        }
        message("\nProcessing lineage: ", lineage_name, " - method: ", method_val)
      } else {
        message("\nProcessing lineage: ", lineage_name, " - method: ", ifelse(use_perpendicular, "perpendicular", "sector"))
      }
      
      # Project to 2D
      result_2d <- project_3d_to_2d(
        sap3d = sap3d_subset,
        lineage_filter = lineage_filter,
        K = K,
        centroids = ap_result$centroids,
        use_perpendicular = current_method,
        ap_sectors = ap_sectors,
        sector_col = sector_col,
        lineage_col = lineage_col
      )
      
      if (!is.null(result_2d)) {
        results[[lineage_name]] <- result_2d
        # Ensure all data frames have same columns
        if (nrow(all_data) == 0) {
          all_data <- result_2d
        } else {
          # Get union of all columns
          all_cols <- unique(c(colnames(all_data), colnames(result_2d)))
          # Ensure both data frames have all columns
          for (col in all_cols) {
            if (!(col %in% colnames(all_data))) {
              all_data[[col]] <- NA
            }
            if (!(col %in% colnames(result_2d))) {
              result_2d[[col]] <- NA
            }
          }
          # Arrange columns in same order
          all_data <- all_data[, all_cols]
          result_2d <- result_2d[, all_cols]
          # Merge data
          all_data <- rbind(all_data, result_2d)
        }
      }
    }
  } else {
    # Use original lineages parameter
    # Step 1: Process all lineages and collect data
    for (lineage_name in names(lineages)) {
      lineage_filter <- lineages[[lineage_name]]
      
      # Determine method for current lineage
      current_method <- use_perpendicular
      if (!is.null(lineage_methods) && lineage_name %in% names(lineage_methods)) {
        method_val <- lineage_methods[[lineage_name]]
        if (method_val == "sector") {
          current_method <- FALSE
        } else if (method_val == "perpendicular") {
          current_method <- TRUE
        }
        message("\nProcessing lineage: ", lineage_name, " (", paste(lineage_filter, collapse = ", "), ") - method: ", method_val)
      } else {
        message("\nProcessing lineage: ", lineage_name, " (", paste(lineage_filter, collapse = ", "), ") - method: ", ifelse(use_perpendicular, "perpendicular", "sector"))
      }
      
      # Project to 2D
      result_2d <- project_3d_to_2d(
        sap3d = sap3d,
        lineage_filter = lineage_filter,
        K = K,
        centroids = ap_result$centroids,
        use_perpendicular = current_method,
        ap_sectors = ap_sectors,
        sector_col = sector_col,
        lineage_col = lineage_col
      )
      
      if (!is.null(result_2d)) {
        results[[lineage_name]] <- result_2d
        # Ensure all data frames have same columns
        if (nrow(all_data) == 0) {
          all_data <- result_2d
        } else {
          # Get union of all columns
          all_cols <- unique(c(colnames(all_data), colnames(result_2d)))
          # Ensure both data frames have all columns
          for (col in all_cols) {
            if (!(col %in% colnames(all_data))) {
              all_data[[col]] <- NA
            }
            if (!(col %in% colnames(result_2d))) {
              result_2d[[col]] <- NA
            }
          }
          # Arrange columns in same order
          all_data <- all_data[, all_cols]
          result_2d <- result_2d[, all_cols]
          # Merge data
          all_data <- rbind(all_data, result_2d)
        }
      }
    }
  }
  
  # Calculate unified coordinate range (if data exists)
  x_range <- NULL
  y_range <- NULL
  if (nrow(all_data) > 0) {
    x_min <- min(all_data$x2d, na.rm = TRUE)
    x_max <- max(all_data$x2d, na.rm = TRUE)
    y_min <- min(all_data$y2d, na.rm = TRUE)
    y_max <- max(all_data$y2d, na.rm = TRUE)
    
    # Add some margin
    x_pad <- (x_max - x_min) * 0.1
    y_pad <- (y_max - y_min) * 0.1
    
    x_range <- c(x_min - x_pad, x_max + x_pad)
    y_range <- c(y_min - y_pad, y_max + y_pad)
    
    message("\nUnified coordinate range:")
    message("  X: [", round(x_range[1], 2), ", ", round(x_range[2], 2), "]")
    message("  Y: [", round(y_range[1], 2), ", ", round(y_range[2], 2), "]")
  }
  
  # Step 2: Save all results (using unified coordinate range)
  for (lineage_name in names(results)) {
    result_2d <- results[[lineage_name]]
    # Save results
    save_results(
      data = result_2d,
      fn = fn0,
      lineage_name = lineage_name,
      save_html = save_html,
      html_template = as.data.frame(ap_result$centroids[, c("x", "y")]),
      x_range = x_range,
      y_range = y_range,
      centroids = ap_result$centroids,
      ap_sectors = ap_sectors
    )
  }
  
  message("\n=== Processing complete ===")
  setwd(original_wd)
  message("Output files saved in: ", normalizePath(output_dir))
  
  # Switch back to original working directory
  message("Switched back to original working directory: ", original_wd)
  
  return(list(
    K = K,
    centroids = ap_result$centroids,
    results = results
  ))
}
