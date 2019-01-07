#'@title [HSU], [HSD], [AHSU] and [AHSD] procedures
#'
#'@description
#'Apply the [HSU], [HSD], [AHSU] and [AHSD] procedures,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'\code{DBH} and \code{ADBH} are wrapper functions for \code{discrete.BH}. 
#'\code{DBH} simply passes all its parameters to \code{discrete.BH} with \code{adaptive = FALSE}. 
#'\code{ADBH} does the same with \code{adaptive = TRUE}.
#'
#'This version: 2018-11-13.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{DiscreteFDR}}, \code{\link{DBR}}
#'
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv FALSE
#'@templateVar pCDFlist TRUE
#'@templateVar raw.pvalues TRUE
#'@templateVar direction TRUE
#'@templateVar ret.crit.consts TRUE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@templateVar adaptive TRUE
#'@template param 
#'
#'@template example
#'@examples
#'
#'DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#'DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#'DBH.sd.fast$Adjusted
#'DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'DBH.sd.crit$Adjusted
#'
#'ADBH.su.fast <- ADBH(raw.pvalues, pCDFlist)
#'ADBH.sd.fast <- ADBH(raw.pvalues, pCDFlist, direction = "sd")
#'ADBH.sd.fast$Adjusted
#'ADBH.su.crit <- ADBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'ADBH.sd.crit <- ADBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'ADBH.sd.crit$Adjusted
#'
#'@templateVar DBR FALSE
#'@template return
#'
#'@name discrete.BH
NULL

#'@rdname discrete.BH
#'@export
discrete.BH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", adaptive = FALSE, ret.crit.consts = FALSE){
  m <- as.integer(length(raw.pvalues))
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
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(direction == "su"){
    # SU case
    # apply the shortcut drawn from Lemma 2, that is
    # c.m >= the effective critical value associated to alpha/(1 + alpha)
    pv.list.c.m <- short.eff(pv.list.all, alpha/(1 + alpha))
    # compute transformed support
    y.c.m <- kernel.DBH.fast(stepf, pv.list.c.m, pv.list.c.m)
    # search the values of the vector <= m*alpha
    Ceiling <- which(y.c.m <= m * alpha)
    # keep the greatest value (note that Ceiling is sorted in nondecreasing order)
    c.m <- pv.list.c.m[Ceiling[length(Ceiling)]]
    
    if(ret.crit.consts){
      # restrict attention to these values, because c.k needs to be <= c.m
      pv.list <- pv.list.all[pv.list.all <= c.m] 
      if(adaptive){
        # apply the shortcut drawn from Lemma 4, that is
        # c.1 >= the effective critical value associated to min((1-c.m)*alpha/m,c.m)
        pv.list <- short.eff(pv.list, min((1 - c.m) * alpha / m, c.m))
        # compute transformed support
        y <- kernel.ADBH.crit(stepf, pv.list, c.m, alpha, NULL) 
      }
      else{
        # apply the shortcut drawn from Lemma 2, that is
        # c.1 >= the effective critical value associated to (alpha/m)/(1 + alpha)
        pv.list <- short.eff(pv.list, alpha/(m + m * alpha))
        # compute transformed support
        y <- kernel.DBH.crit(stepf, pv.list, c.m, alpha, NULL)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      # restrict attention to these values, because rejected p-values need to be <= c.m
      obs.pvals <- sorted.pvals[sorted.pvals <= c.m]
      if(length(obs.pvals)){
        if(adaptive){
          # compute transformed observed p-values
          y <- kernel.ADBH.fast(stepf, obs.pvals, c.m)
        }
        else{
          # compute transformed observed p-values
          y <- kernel.DBH.fast(stepf, obs.pvals, c.m)
        }
        # determine significant (transformed) p-values
        idx <- which(y <= 1:length(obs.pvals) * alpha)
      }
      else{
        idx <- integer(0)
      }
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      m.rej <- 0
      idx <- integer(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    # SD case
    if(ret.crit.consts){
      # apply the shortcut drawn from Lemma 3, that is
      # c.1 >= the effective critical value associated to (alpha/m)/(1 + alpha/m)
      pv.list <- short.eff(pv.list.all, alpha/(m + alpha))
      # then re-add the observed p-values (needed to compute the adjusted p-values),
      # because we may have removed some of them the shortcut
      pv.list <- sort(unique(c(pv.list, sorted.pvals)))
      if(adaptive){
        # compute transformed support
        y <- kernel.ADBH.crit(stepf, pv.list, pv.list, alpha, sorted.pvals)
      }
      else{
        # compute transformed support
        y <- kernel.DBH.crit(stepf, pv.list, pv.list, alpha, sorted.pvals)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      if(adaptive){
        # compute transformed sorted p-values
        y <- kernel.ADBH.fast(stepf, sorted.pvals, sorted.pvals)
      }
      else{
        # compute transformed sorted p-values
        y <- kernel.DBH.fast(stepf, sorted.pvals, sorted.pvals)
      }
      # determine significant (transformed) p-values
      idx <- which(y > 1:m * alpha) 
    }
    if(length(idx)){
      m.rej <- min(idx) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- raw.pvalues[idx]
      }
      else{
        idx <- numeric(0)
        pvec.rej <- numeric(0)
      }
    }
    else{
      m.rej <- m
      idx <- 1:m
      pvec.rej <- raw.pvalues
    }
  }
  #--------------------------------------------
  #       Create output list
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, k.hat = m.rej)
  if(direction == "sd"){
    if(ret.crit.consts){
      # compute adjusted p-values
      pv.adj = cummax(pmin(y$pval.transf / 1:m, 1))
    }
    else{
      # compute adjusted p-values
      pv.adj = cummax(pmin(y / 1:m, 1))
    }
    # add adjusted p-values to output list
    ro <- order(o)
    output <- c(output, list(Adjusted = pv.adj[ro]))  
  }
  # add critical constants to output list
  if(ret.crit.consts) output$Critical.constants = crit.constants
  return(output)
}

#'@rdname discrete.BH
#'@export
DBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = FALSE, ret.crit.consts))
}

#'@rdname discrete.BH
#'@export
ADBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = TRUE, ret.crit.consts))
}