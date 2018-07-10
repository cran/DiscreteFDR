#'@title Kernel functions for DBH 
#'
#'@description
#'Kernel functions that transform 
#'observed p-values or their support according to [HSU] and [HSD].
#'The output is used by \code{\link{DBH}}.
#'The end user should not use it.
#'
#'@details
#'This version: 2017-09-14.
#'
#'@seealso
#'\code{\link{DBH}}, \code{\link{DiscreteFDR}}
#'
#'@templateVar msg TRUE
#'@templateVar stepf TRUE
#'@templateVar pv.numer TRUE
#'@templateVar pv.denom TRUE
#'@templateVar alpha FALSE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem TRUE
#'@templateVar verbose TRUE
#'@templateVar pCDFlist FALSE
#'@templateVar raw.pvalues FALSE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@template param 
#'
#'@template example
#'@examples
#'
#'m <- length(raw.pvalues)
#'alpha <- 0.05
#'
#'#Compute the step functions from the supports
#'stepf <- build.stepfuns(pCDFlist)
#'
#'#We stay in a step-down context, where pv.numer=pv.denom,
#'#for the sake of simplicity
#'
#'#If not searching for critical constants, we use only the observed p-values
#'sorted.pvals <- sort(raw.pvalues)
#'y <- kernel.DBH("", stepf, sorted.pvals, sorted.pvals)
#'
#'#If searching for critical constants, we use (almost) the complete support
#'pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
#'# apply the shortcut drawn from Lemma 3, that is
#'# c.1 >= the effective critical value associated to (alpha/m)/(1 + alpha/m)
#'pv.list<-short.eff(pv.list.all, (alpha/m)/(1 + alpha/m) )
#'# then re-add the observed p-values (needed to compute the adjusted p-values),
#'# because we may have removed some of them the shortcut
#'pv.list <- sort(unique(c(pv.list, sorted.pvals)))
#'# compute transformed support
#'y <- kernel.DBH("", stepf, pv.list, pv.list)
#'
#'@return
#'A vector of transformed p-values.
#'
#'@name kernel.DBH
#'@export
kernel.DBH <- function(msg = "", stepf, pv.numer, pv.denom, bigMem = FALSE, verbose = TRUE){
  m <- length(stepf)
  # the output is the vector y= \sum_{i=1}^m F_i(pv.numer)/(1-F_i(pv.denom))
  # note that pv.number is either all values of all supports of the p-values
  # (when searching for the critical values), or only the observed p-values
  # (when skipping the critical values computation)
  if(bigMem){
    # the small data size allows to use lapply
    if(verbose) message(paste(msg, "Using Reduce ...", sep = ""))
    if(all(pv.denom == pv.numer)){
      # SD case, see (11)
      y <- tryCatch(Reduce('+', lapply(stepf, function(s){v <- s(pv.numer);return(v / (1 - v))})),
                    error = function(err){return(kernel.DBH(err, stepf, pv.numer, pv.numer))})
    }else{
      # SU case, see (10)
      y <- tryCatch(Reduce('+', lapply(stepf, function(s){return(s(pv.numer) / (1 - s(pv.denom)))})),
                    error = function(err){return(kernel.DBH(err, stepf, pv.numer, pv.denom))})
    }
  }
  else{
    # the large data size requires to use a for loop
    if(verbose && msg != "") message(paste(msg, "Using loop ...", sep = ""))
    y <- 0
    if(all(pv.denom == pv.numer)){
      # SD case, see (11)
      for(i in 1:m){
        s <- stepf[[i]]
        s_numer <- s(pv.numer)
        y <- y + s_numer / (1 - s_numer)
      }
    }else{
      for(i in 1:m){
        # SU case, see (10)
        s <- stepf[[i]]
        s_numer <- s(pv.numer)
        s_denom <- 1 - s(pv.denom)
        y <- y + s_numer / s_denom
      }
    }
  }
  return(y)
}

#'@title [HSU] and [HSD] procedures
#'
#'@description
#'Apply the [HSU] or [HSD] procedure,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'This version: 2018-02-20.
#'
#'@templateVar msg FALSE
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem TRUE
#'@templateVar verbose TRUE
#'@templateVar pCDFlist TRUE
#'@templateVar raw.pvalues TRUE
#'@templateVar direction TRUE
#'@templateVar ret.crit.consts TRUE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@template param 
#'
#'@template example
#'@examples
#'
#'DBH.su.fast <- DBH(pCDFlist, raw.pvalues)
#'DBH.sd.fast <- DBH(pCDFlist, raw.pvalues, direction="sd")
#'DBH.su.crit <- DBH(pCDFlist, raw.pvalues, ret.crit.consts=TRUE)
#'DBH.sd.crit <- DBH(pCDFlist, raw.pvalues, direction="sd", ret.crit.consts=TRUE)
#'
#'@templateVar DBR FALSE
#'@template return
#'
#'@name DBH
#'@export
DBH <- function(pCDFlist, raw.pvalues, alpha = 0.05, direction = "su", ret.crit.consts = FALSE, bigMem = FALSE, verbose = FALSE){
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
  pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(direction == "su"){
    # SU case
    # apply the shortcut drawn from Lemma 2, that is
    # c.m >= the effective critical value associated to alpha/(1 + alpha)
    pv.list.c.m<-short.eff(pv.list.all, alpha/(1 + alpha) )
    # compute transformed support
    y.c.m <- kernel.DBH("", stepf, pv.list.c.m, pv.list.c.m, bigMem, verbose)
    # search the values of the vector <= m*alpha
    Ceiling<-which(y.c.m<=m*alpha)
    # keep the greatest value (note that Ceiling is sorted in nondecreasing order)
    c.m<-pv.list.c.m[Ceiling[length(Ceiling)]]
    rm(pv.list.c.m,y.c.m,Ceiling)
    if(ret.crit.consts){
      # apply the shortcut drawn from Lemma 2, that is
      # c.1 >= the effective critical value associated to (alpha/m)/(1 + alpha)
      pv.list<-short.eff(pv.list.all, (alpha/m)/(1 + alpha) )
      # restrict attention to these values, because c.k needs to be <= c.m
      pv.list <- pv.list[pv.list <= c.m] 
      # compute transformed support
      y <- kernel.DBH("", stepf, pv.list, c.m, bigMem, verbose)
      # find critical constants
      crit.constants <- c(sapply(1:(m - 1), function(i){pv.list[which.max(y[y <= i * alpha])]}), c.m)
      idx <- which(sorted.pvals <= crit.constants)
    }else{
      # restrict attention to these values, because rejected p-values need to be <= c.m
      obs.pvals<-sorted.pvals[sorted.pvals<=c.m]
      # compute transformed observed p-values
      y <- kernel.DBH("", stepf, obs.pvals, c.m, bigMem, verbose)
      # determine significant (transformed) p-values
      idx <- which(y <= seq_len(length(obs.pvals)) * alpha)
      rm(obs.pvals)
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[idx]
    }else{
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
      pv.list<-short.eff(pv.list.all, (alpha/m)/(1 + alpha/m) )
      # then re-add the observed p-values (needed to compute the adjusted p-values),
      # because we may have removed some of them the shortcut
      pv.list <- sort(unique(c(pv.list, sorted.pvals)))
      # compute transformed support
      y <- kernel.DBH("", stepf, pv.list, pv.list, bigMem, verbose)
      # find critical constants
      crit.constants <- sapply(1:m, function(i){pv.list[which.max(y[y <= i * alpha])]})
      idx <- which(sorted.pvals > crit.constants)
    }else{
      # compute transformed sorted p-values
      y <- kernel.DBH("", stepf, sorted.pvals, sorted.pvals, bigMem, verbose)
      # determine significant (transformed) p-values
      idx <- which(y > 1:m * alpha) 
    }
    if(length(idx)){
      m.rej <- min(idx) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- raw.pvalues[idx]
      }else{
        idx <- numeric(0)
        pvec.rej <- numeric(0)
      }
    }else{
      m.rej <- m
      idx <- 1:m
      pvec.rej <- raw.pvalues
    }
  }
  #--------------------------------------------
  #       Create output list
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, Max.k = m.rej)
  if(direction == "sd"){
    if(ret.crit.consts){
      # get the indexes of the observed p-values inside pv.list
      idx <- sapply(1:m, function(i){which(pv.list == sorted.pvals[i])})
      # compute adjusted p-values
      pv.adj = cummax(pmin(y[idx] / 1:m, 1))
    }
    else{
      # compute adjusted p-values
      pv.adj = cummax(pmin(y / 1:m, 1))
    }
    # add adjusted p-values to output list
    ro <- order(o)
    output <- c(output, list(Adjusted = pv.adj[ro]))  
  }
  if(ret.crit.consts) output <- c(output, list(Critical.constants = crit.constants)) # add critical constants to output list
  return(output)
}