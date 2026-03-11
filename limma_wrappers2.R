library(dplyr)


#' Differential Expression Wrapper for Limma
#'
#' A wrapper function that performs differential expression analysis using the \code{limma} 
#' pipeline (lmFit and eBayes). It processes the results into a tidy format, calculates 
#' descriptive statistics, and identifies significant hits based on a specified FDR threshold.
#'
#' @param design_mat A design matrix created by \code{model.matrix} representing the 
#' experimental design.
#' @param expr_data A numeric matrix or data frame containing expression values, 
#' with genes in rows and samples in columns.
#' @param coefs An integer or vector specifying which coefficients (columns) of the 
#' design matrix to test for differential expression. Default is 2.
#' @param threshold A numeric value representing the False Discovery Rate (FDR) 
#' cutoff for significance. Default is 0.01.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{FULL_RESULTS}: A tibble containing the topTable output, additional 
#'   gene IDs, significance flags, and descriptive statistics.
#'   \item \code{HITS}: A character vector of gene names that passed the FDR threshold.
#'   \item \code{FIT}: The \code{MArrayLM} object returned by \code{eBayes}.
#' }
#' 
#' @import limma
#' @import dplyr
#' @import tibble
#' @import pastecs

limma_wrapp = function(design_mat, expr_data, coefs = c(2), threshold = 0.01) {
  
  # Fit Data
  diffExpr <- limma::lmFit(expr_data, design_mat)
  fit.cont <- limma::eBayes(diffExpr)
  
  # Get Summary of Pvalues and logFCs
  topDiff <- limma::topTable(fit.cont, coef = coefs, adjust = "fdr", number = Inf)
  
  message(paste("TOTAL HITS:", sum(topDiff$adj.P.Val < threshold)))
  
  # Format for export
  topDiff <- tibble::rownames_to_column(topDiff, "Gene")
  
  # Note: extract_geneID must be defined elsewhere in your package/environment
  topDiff$Gene2 <- extract_geneID(topDiff$Gene)
  
  names(topDiff)[6] <- 'FDR'
  topDiff$Is.FDR.Significant <- ifelse(topDiff$FDR < threshold, 1, 0)
  
  # Move Gene2 column and calculate rank
  topDiff <- topDiff %>% 
    dplyr::relocate(Gene2, .before = logFC) %>%
    dplyr::mutate(PCT_RANK_FDR = dplyr::percent_rank(FDR))
  
  # Calculate descriptive statistics
  desc_ <- t(pastecs::stat.desc(t(expr_data)))
  desc <- desc_[topDiff$Gene, ]
  
  # Identify significant genes
  limmaHits <- topDiff$Gene[topDiff$FDR < threshold]
  
  return(list(
    `FULL_RESULTS` = topDiff, 
    `HITS` = limmaHits,
    `FIT` = fit.cont
  ))
}