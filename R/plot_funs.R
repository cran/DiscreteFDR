#'@name hist.DiscreteFDR
#'@title Histogram of Raw p-Values
#'
#'@description
#'Computes a histogram of the raw p-values of a \code{DiscreteFDR} object.
#'
#'@param x          an object of class "\code{DiscreteFDR}".
#'@param breaks     as in \code{\link{hist}}; here, the Friedman-Diaconis
#'                  algorithm(\code{"FD"}) is used as default.
#'@param ...        further arguments to \code{\link{hist}} or
#'                  \code{\link{plot.histogram}}, respectively.
#'
#'@details
#'This method simply calls \code{\link{hist}} and passes the raw p-values of
#'the object.
#'
#'@return
#'An object of class \code{histogram}.
#'
#'@template example
#'@examples
#'
#'DBH <- DBH(raw.pvalues, pCDFlist)
#'hist(DBH)
#'
#'@importFrom graphics hist
#'@export
hist.DiscreteFDR <- function(x, breaks = "FD", ...){
  # necessary to appease automated R CMD check on CRAN
  main <- xlab <- NULL
  
  # Main title
  if(hasArg(main)) titl <- main else titl <- "Histogram of raw p-values"
  
  # X-axis label
  if(hasArg(xlab)) labl <- xlab else labl <- "Raw p-values"
  
  # call 'hist' function with raw p.values (default breaks: "FD"); "..." passes
  # all additional 'hist' arguments to this call
  hist(x$Data$raw.pvalues, breaks = breaks, main = titl, xlab = labl,  ...)
}


#'@name plot.DiscreteFDR
#'@title Plot Method for \code{DiscreteFDR} objects
#'
#'@description
#'Plots raw p-values of a \code{DiscreteFDR} object and highlights rejected and
#'accepted p-values. If present, the critical values are plotted, too.
#'
#'@param x          an object of class "\code{DiscreteFDR}".
#'@param col        a numeric or character vector of length 3 indicating the
#'                  colors of the \enumerate{
#'                    \item rejected p-values
#'                    \item accepted p-values
#'                    \item critical values (if present).
#'                  }
#'@param pch        a numeric or character vector of length 3 indicating the
#'                  point characters of the \enumerate{
#'                    \item rejected p-values
#'                    \item accepted p-values
#'                    \item critical values (if present and \code{type.crit}
#'                          is a plot type like \code{'p'}, \code{'b'} etc.).
#'                  }
#'@param lwd        a numeric vector of length 3 indicating the thickness of
#'                  the points and lines.
#'@param type.crit  1-character string giving the type of plot desired for the
#'                  critical values (e.g.: \code{'p'}, \code{'l'} etc; see
#'                  \code{\link{plot}}).
#'@param legend     if NULL, no legend is plotted; otherwise expecting a
#'                  character string like "topleft" etc. or a numeric vector
#'                  of two elements indicating (x, y) coordinates.
#'@param ...        further arguments to \code{\link{plot.default}}.
#'
#'@template example
#'@examples
#'
#'DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#'DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#'DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'
#'plot(DBH.sd.fast)
#'plot(DBH.sd.crit, xlim = c(1, 5), ylim = c(0, 0.4))
#'plot(DBH.su.fast, col = c(2, 4), pch = c(2, 3), lwd = c(2, 2), 
#'     legend = "topleft", xlim = c(1, 5), ylim = c(0, 0.4))
#'plot(DBH.su.crit, col = c(2, 4, 1), pch = c(1, 1, 4), lwd = c(1, 1, 2), 
#'     type.crit = 'o', legend = c(1, 0.4), lty = 1, xlim = c(1, 5), 
#'     ylim = c(0, 0.4))
#'
#'@importFrom graphics plot lines points
#'@importFrom methods hasArg
#'@export
plot.DiscreteFDR <- function(x, col = c(2, 4, 1), pch = c(1, 1, 1), lwd = c(1, 1, 1), type.crit = 'b', legend = NULL, ...){
  # necessary to appease automated R CMD check on CRAN
  main <- ylab <- lty <- NULL
  
  # determine number of tests and rejections
  m <- length(x$Data$raw.pvalues)
  k <- length(x$Indices)
  
  # replicate shorter plot parameter vectors to avoid errors
  col <- rep_len(col, 3)
  pch <- rep_len(pch, 3)
  lwd <- rep_len(lwd, 3)
  
  # Main title
  if(hasArg(main)) titl <- main else titl <- x$Method
  
  # Y-axis label
  if(hasArg(ylab)) labl <- ylab else labl <- "Sorted raw p-values"
  
  # start plotting with empty area
  plot(x$Data$raw.pvalues, col = NA, main = titl, ylab = labl, ...)
  # plot rejected p-values
  points(sort(x$Data$raw.pvalues)[1:k], col = col[1], pch = pch[1], lwd = lwd[1], ...)
  # plot accepted p-values
  points((k + 1):m, sort(x$Data$raw.pvalues[-x$Indices]), col = col[2], pch = pch[2], lwd = lwd[2], ...)
  
  # plot critical values (if present and plotting is requested by the user)
  if(!is.null(x$Critical.values) && type.crit != 'n'){
    lines(x$Critical.values, col = col[3], lwd = lwd[3], pch = pch[3], type = type.crit, ...)
  }
  
  # plot legend
  if(!is.null(legend)){
    n <- 3 - is.null(x$Critical.values)
    lt <- rep(0, n)
    if(n > 2 && hasArg(lty) && !(type.crit %in% c('p', 'n'))) lt[3] <- list(...)$lty
    if(length(legend) == 1){
      legend(legend, NULL, c("Rejected", "Accepted", "Critical values")[1:n], col = col[1:n], pch = pch[1:n], lty = lt[1:n], lwd = lwd[1:n])
    }else if(length(legend) == 2){
      legend(legend[1], legend[2], c("Rejected", "Accepted", "Critical values")[1:n], col = col[1:n], pch = pch[1:n], lty = lt[1:n], lwd = lwd[1:n])
    }else warning("Expecting character string or numeric vector of one or two elements for creating a legend")
  }
}
