#' Variance-based row filtering
#'
#' @export
rna.mafilter <- function(expr, rownames = NULL, var.cutoff = 0.5)
{
    vars <- apply(expr, 1, IQR)
    if (0 < var.cutoff && var.cutoff < 1) {
      quant <- quantile(vars, probs = var.cutoff)
      var.keep <- !is.na(vars) & vars > quant
      res <- expr[var.keep,]
    }
    else stop("Cutoff Quantile has to be between 0 and 1.")


    if(!is.null(rownames)){
      na.keep <- !is.na(rownames)

      res <- expr[na.keep,]
      rownames <- rownames[na.keep]

      res <- aggregate(res, list(rownames), mean)
      rownames(res) <- res[,"Group.1"]
      res <- res[,-1]
    }

    return(res)
}
