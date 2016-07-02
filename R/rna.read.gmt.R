#' Read a GMT file
#'
#' Read a gene set (.gmt) file into R list
#'
#' @param file  text file
#'
#' @export

rna.read.gmt <- function(file){
    f <- readLines(file)
    lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
    names(lst) = sapply(lst, function(x) x[1])
    lst = lapply(lst, function(x) x[-(1:2)])
    return(lst)
}
