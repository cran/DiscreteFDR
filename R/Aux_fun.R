#' @name match.pvals
#' 
#' @title
#' Matching Raw P-Values with Supports 
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Constructs the observed p-values from the raw observed p-values, by rounding
#' them to their nearest neighbor matching with the supports of their
#' respective CDFs (as in function `p.discrete.adjust()` of package
#' `discreteMTP`, which is no longer available on CRAN).
#' 
#' **Note**: In the next version, this is to become an internal function and
#' will have to be called directly via `:::`, i.e.
#' `DiscreteFDR:::match.pvals()`.
#' 
#' @details
#' Well computed raw p-values should already belong to their respective CDF
#' support. So this function is called at the beginning of [discrete.BH()],
#' [DBH()], [ADBH()] and [DBR()], just in case raw p-values are biased.
#'
#' For each raw p-value that needs to be rounded, a warning is issued.
#'
#' @seealso
#' [`discrete.BH()`], [`DBR()`]
#'
#' @templateVar pCDFlist TRUE
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist.indices TRUE
#' @template param
#' 
#' @examples \dontrun{
#' toyList <- list(c(0.3,0.7,1),c(0.1,0.65,1))
#' toyRaw1 <- c(0.3,0.65)
#' match.pvals(toyList,toyRaw1)
#' toyRaw2 <- c(0.31,0.6)
#' match.pvals(toyList,toyRaw2)
#' }
#'
#' @return
#' A vector where each raw p-value has been replaced by its nearest neighbor, if
#' necessary.
#'
#' @importFrom lifecycle deprecate_soft
#' @export
match.pvals <- function(pCDFlist, raw.pvalues, pCDFlist.indices = NULL){
  deprecate_soft("2.0.0", "match.pvals()", 
                 details = paste0("This will become a purely internal",
                                  "function. Please call it directly via",
                                  "'DiscreteFDR:::match.pvals()'.")
                 )
  
  m <- length(raw.pvalues)
  if(!is.null(pCDFlist.indices)){
    idx <- unlist(pCDFlist.indices)
    counts <- sapply(pCDFlist.indices, length)
    pCDFlist <- rep(pCDFlist, counts)[order(idx)]
  }
  n <- length(pCDFlist)
  if(m > 0 && m == n){
    pvec <- raw.pvalues
    in.CDF <- numeric(m)
    for (k in 1:m) {
      in.CDF[k] <- match(pvec[k], pCDFlist[[k]])
      if (is.na(in.CDF[k])){
        in.CDF[k] <- which.min(abs(pCDFlist[[k]] - pvec[k]))
        pvec[k] <- pCDFlist[[k]][in.CDF[k]]
        ordinal <- "th"
        if (k%%10==1) ordinal <- "st"
        if (k%%10==2) ordinal <- "nd"
        if (k%%10==3) ordinal <- "rd"
        if (k%%100-k%%10==10) ordinal <- "th"
        warning("Since ", raw.pvalues[k], 
                " is not a value of the CDF of the ",k,ordinal ," p-value,\n  the p-value is rounded to be ", 
                pCDFlist[[k]][in.CDF[k]], call. = F)
      }
    }
    return(pvec)
  }else{
    stop("'pCDFlist' and 'raw.pvalues' must have the same non-zero length")
  }
}
