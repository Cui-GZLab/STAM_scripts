library(dplyr)

cal_APorder <- function(meta) {
  meta$section <- as.numeric(as.character(meta$section))
  
  A_meta <- filter(meta, AP == "A")
  
  if (nrow(A_meta) == 0) {
    warning("No 'A' region found in the current group.")
    meta$APorder <- NA
    return(meta)
  }
  
  min_A <- min(A_meta$section, na.rm = TRUE)
  max_A <- max(A_meta$section, na.rm = TRUE)
  
  # A ：abs(section - max_A) + 1
  # P ：(section - min_A) + abs(min_A - max_A) + 2
  meta$APorder <- ifelse(meta$AP == "A",
                         abs(meta$section - max_A) + 1,
                         (meta$section - min_A) + abs(min_A - max_A) + 2)
  
  return(meta)
}

