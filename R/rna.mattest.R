#' Array-based t-test
#'
#' Perform two-sample t-test to evaluate differential expression of genes from
#' two experimental conditions or phenotypes
#'
#' @param expr    numeric matrix or data.frame consisting of gene expression values
#' @param ref     numeric vector. Index of reference samples in the array
#' @param samp    numeric vector. Index of case samples in the array
#' @param p.value double or NULL (default). If NULL returns a vector of p-values
#'                for each genes, otherwise returns a logical vecotr of genes with
#'                p-values below the specified threshold.
#' @export

rna.mattest <- function(expr, ref, samp, p.value = NULL){

  ttest <- apply(expr, 1, function(x) t.test(x[ref], x[samp])$p.value)

  if(is.null(p.value))
    return(ttest)

  return(ttest <= p.value)
}
