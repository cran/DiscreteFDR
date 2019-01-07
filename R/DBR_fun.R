#'@title [DBR-lambda] procedure
#'
#'@description
#'Apply the [DBR-lambda] procedure,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'
#'[DBR-lambda] is the discrete version of the [Blanchard-Roquain-lambda] procedure (see References),
#'the authors of the latter suggest to take \code{lambda = alpha} (see their Proposition 17),
#'which explains the choice of the default value here. 
#'
#'This version: 2018-11-13.
#'
#'@section References:
#'G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control under independence and dependence. Journal of Machine Learning Research, 10, 2837-2871.
#'
#'@seealso
#'\code{\link{discrete.BH}}, \code{\link{DiscreteFDR}}
#'
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv FALSE
#'@templateVar pCDFlist TRUE
#'@templateVar raw.pvalues TRUE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts TRUE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda TRUE
#'@template param 
#'
#'@template example
#'@examples
#'
#'DBR.su.fast <- DBR(raw.pvalues, pCDFlist)
#'DBR.su.crit <- DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'
#'@templateVar DBR TRUE
#'@template return
#'
#'@name DBR
#'@export
DBR <- function(raw.pvalues, pCDFlist, alpha = 0.05, lambda = NULL, ret.crit.consts = FALSE){
  m <- length(raw.pvalues)
  # if lambda is not provided (lambda = NULL),
  # set lambda = alpha
  if (is.null(lambda)){
    lambda <- alpha
  }
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(pCDFlist, raw.pvalues)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)                                                                                       
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       define step functions corresponding to distribution functions of p-values under Null-Hypotheses
  #--------------------------------------------
  stepf <- build.stepfuns(pCDFlist)
  #--------------------------------------------
  pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [DBR-lambda] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(ret.crit.consts){
    # construct the vector of all values of all supports of the p-values
    pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
    # apply the shortcut drawn from Corollary 3, that is
    # c.1 >= the effective critical value associated to min((1 - lambda) * alpha/m , lambda)
    pv.list <- short.eff(pv.list.all, min((1 - lambda) * alpha/m , lambda))
    # then re-add the observed p-values (needed to compute the adjusted p-values),
    # because we may have removed some of them the shortcut
    pv.list <- sort(unique(c(pv.list, sorted.pvals)))
    # compute transformed support
    y <- kernel.DBR.crit(stepf, pv.list, lambda, alpha, sorted.pvals)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals <= crit.constants)   
  }
  else{
    # compute transformed p-values
    y <- kernel.DBR.fast(stepf, sorted.pvals, lambda)
    idx <- which( y <= alpha)
  }
  m.rej <- length(idx)
  if(m.rej){
    idx <- which(pvec <= sorted.pvals[m.rej]) 
    pvec.rej <- raw.pvalues[idx]
  }else{
    idx <- integer(0)
    pvec.rej <- numeric(0)
  }
  #--------------------------------------------
  #       Create output list
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, k.hat = m.rej, Lambda = lambda)
  if(ret.crit.consts){
    output <- c(output, list(Critical.constants = crit.constants))
    # compute adjusted p-values
    # recall that transformed p-values where max_i F_i(p) > lambda,
    # that is for indices > y$m.lambda, are set to 1
    pv.adj <- rev(cummin(rev(pmin(y$pval.transf / c(seq_len(y$m.lambda), rep(1, m - y$m.lambda)), 1))))
  }
  else{
    # compute adjusted p-values
    pv.adj <- rev(cummin(rev(pmin(y, 1))))
  }
  # add adjusted p-values to output list
  ro <- order(o)
  output <- c(output, list(Adjusted = pv.adj[ro]))
  return(output)
}