#' Pathway t-test
#'
#' Perform two-sample t-test to evaluate differential expression of pathways
#' obtained from a gene set
#'
#' @param expr        numeric matrix or data.frame consisting of gene expression values
#' @param gs          list. gene set collection
#' @param genelist    character vector. Names of genes in the dataset.
#'                    Lenght of genelist and n rows of expr must match
#' @param ref         numeric vector. Index of reference samples in the array
#' @param samp        numeric vector. Index of case samples in the array
#' @param set.size    vector. First and second elements specify the minimum and
#'                    maximum numbers of gene per pathway respectively
#' @export

rna.pathwayttest<- function(expr, gs, genelist, ref, samp, set.size = c(10,500)){

  gs <- gs[(sapply(gs, length) >= set.size[1]) & (sapply(gs, length) <= set.size[2])]

  gs.means <- lapply(gs, function(x) colMeans(expr[match(x, genelist),], na.rm = TRUE))
  gs.ttest <- lapply(gs.means, function(x){ t <- t.test(x[samp], x[ref])
                                  data.frame(t$p.value, t$statistic)
                                  })
  pathway.names <- names(gs.ttest)
  gs.ttest <- do.call(rbind, gs.ttest)

  res <- data.frame(collection = pathway.names, p.value = gs.ttest[,1], fdr = p.adjust(gs.ttest[,1], "BH"), t.value = gs.ttest[,2], logFC = rna.foldchange(do.call(rbind, gs.means), ref = ref, samp = samp), members = get.members(gs, pathway.names))



  res <- res[order(res$p.value, decreasing = FALSE),]
}

get.members <- function(geneset, pathways){
  ix <- match(pathways, names(geneset))

  glist <- sapply(ix, function(x) toString(paste(geneset[[x]], collapse = "///")))
  return(glist)
}
