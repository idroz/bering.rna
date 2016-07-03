#' Compute log2-fold changes
#'
#' Fold changes are computed between sample and reference arrays
#'
#' @param expr    numeric matrix or data.frame consisting of gene expression values
#' @param ref     numeric vector. Index of reference samples in the array
#' @param samp    numeric vector. Index of case samples in the array
#'
#' @export

rna.foldchange <- function(expr, ref, samp){

  ref.means <- rowMeans(expr[,ref], na.rm = TRUE)
  samp.means <- rowMeans(expr[,samp], na.rm = TRUE)

  fc <- samp.means/ref.means

  # Fold changes < 1 become negative
  logFC <- log2(fc)

  return(sign(logFC) * (2^logFC))

}
