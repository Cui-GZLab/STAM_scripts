#' VGLM-based Differential Expression Analysis
#' 
#' Performs differential expression analysis using VGAM::vglm with natural 
#' spline smoothing to test for expression changes across a continuous factor.
#' 
#' @param seu_obj A Seurat object containing gene expression data and metadata
#' @param factor_name Name of the metadata column to test against. 
#'   Common values include "section" or "APoder" (AP axis). 
#'   (default: "section")
#' @param df Degrees of freedom for natural spline (default: 3)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#' 
#' @return A data frame with columns:
#'   - pval: Raw p-values from vglm test
#'   - qval: Benjamini-Hochberg adjusted p-values
#' 
#' @examples
#' \dontrun{
#' select_seu <- readRDS("E5.5_Epi_select_seu_obj_gene_filt.rds")
#' deg_result <- vglm_smooth_deg(select_seu, factor_name = "section")
#' write.csv(deg_result, "vglm_deg_result.csv")
#' }
#' 
#' @export
vglm_smooth_deg <- function(seu_obj, factor_name = "section", df = 3, verbose = FALSE) {
  
  count_matrix <- seu_obj@assays$RNA@counts
  genes <- rownames(count_matrix)
  factor_dat <- seu_obj@meta.data[, factor_name]
  
  test_gene <- function(gene) {
    data <- data.frame(exp = log10(count_matrix[gene, ] + 1),
                       factor = factor_dat)
    colnames(data) <- c("exp", factor_name)
    
    full_model_fit <- VGAM::vglm(
      as.formula(paste0("exp ~ sm.ns(", factor_name, ", df=", df, ")")),
      data = data,
      epsilon = 1e-1,
      family = uninormal()
    )
    
    pval <- coef(summary(full_model_fit))[3, 4]
    return(pval)
  }
  
  if (verbose) {
    cat(paste("Testing", length(genes), "genes against", factor_name, "\n"))
  }
  
  pval <- sapply(genes, test_gene)
  qval <- p.adjust(pval, method = "BH")
  
  result <- data.frame(pval = pval, qval = qval)
  return(result)
}
