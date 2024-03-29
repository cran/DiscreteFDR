#' @title
#' Matching raw p-values with supports 
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Constructs the observed p-values from the raw observed p-values, by rounding
#' them to their nearest neighbor matching with the supports of their
#' respective CDFs (as in function `p.discrete.adjust` of package `discreteMTP`,
#' which is no longer available on CRAN).
#' The end user should not use it directly.
#' 
#' **Note**: In future versions, this function will no longer be exported to the
#' global namespace. Instead, it will be a purely internal function and will
#' have to be called directly via `:::`, i.e. `DiscreteFDR:::match.pvals()`.
#' 
#' @details
#' Well computed raw p-values should already belong to their respective CDF
#' support. So this function is called at the beginning of [discrete.BH],
#' [DBH], [ADBH] and [DBR], just in case raw p-values are biased.
#'
#' For each raw p-value that needs to be rounded, a warning is issued.
#'
#' @seealso
#' [discrete.BH], [DBR]
#'
#' @templateVar stepf FALSE
#' @templateVar pv.numer FALSE
#' @templateVar pv.denom FALSE
#' @templateVar alpha FALSE
#' @templateVar sorted.pv FALSE
#' @templateVar pCDFlist TRUE
#' @templateVar raw.pvalues TRUE
#' @templateVar direction FALSE
#' @templateVar ret.crit.consts FALSE
#' @templateVar lambda FALSE
#' @template param
#' 
#' @examples
#' toyList <- list(c(0.3,0.7,1),c(0.1,0.65,1))
#' toyRaw1 <- c(0.3,0.65)
#' match.pvals(toyList,toyRaw1)
#' toyRaw2 <- c(0.31,0.6)
#' match.pvals(toyList,toyRaw2)
#'
#' @return
#' A vector where each raw p-value has been
#' replaced by its nearest neighbor.
#'
#' @name match.pvals
#' @importFrom lifecycle deprecate_soft
#' @export
match.pvals <- function(pCDFlist, raw.pvalues){
  deprecate_soft("1.3.7", "match.pvals()",
                 details = paste("In future versions of this package, this",
                                 "function will no longer be exported to the",
                                 "global namespace and will have to be",
                                 "accessed directly via ':::', i.e.",
                                 "'DiscreteFDR:::match.pvals()'."))
  
  m <- length(raw.pvalues)
  if(m > 0){
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
    return(raw.pvalues)
  }
}
