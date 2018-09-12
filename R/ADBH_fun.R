#'@title Kernel functions for ADBH 
#'
#'@description
#'Kernel functions that transform 
#'observed p-values or their support according to [AHSU] and [AHSD].
#'The output is used by \code{\link{ADBH}}.
#'Additionally, \code{kernel.ADBH.crit} computes and returns
#'the critical constants.
#'The end user should not use them.
#'
#'@details
#'When computing critical constants under step-down,
#'that is, when using \code{kernel.ADBH.crit} with
#'\code{pv.numer=pv.denom}, we still need
#'to get transformed p-values to compute
#'the adjusted p-values. Also, note that here the critical constants are computed by
#'the kernel function and not by the principal
#'function \code{\link{ADBH}}, contrary to what happens
#'with \code{\link{DBH}}. This is why \code{sorted.pv} is needed.
#'
#'This version: 2018-02-20.
#'
#'@seealso
#'\code{\link{ADBH}}, \code{\link{DiscreteFDR}},
#'\code{\link{kernel.DBH}}, \code{\link{DBH}}
#'
#'@templateVar msg TRUE
#'@templateVar stepf TRUE
#'@templateVar pv.numer TRUE
#'@templateVar pv.denom TRUE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv TRUE
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
#'y <- kernel.ADBH.fast("", stepf, sorted.pvals, sorted.pvals)
#'
#'#If searching for critical constants, we use (almost) the complete support
#'pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
#'# apply the shortcut drawn from Lemma 4, that is
#'# c.1 >= the effective critical value associated to alpha/(m + alpha)
#'pv.list<-short.eff(pv.list.all, alpha/(m + alpha) )
#'# then re-add the observed p-values (needed to compute the adjusted p-values),
#'# because we may have removed some of them the shortcut
#'pv.list <- sort(unique(c(pv.list, sorted.pvals)))
#'# compute transformed support
#'y <- kernel.ADBH.crit("", stepf, pv.list, pv.list, alpha, sorted.pvals)
#'crit.constants <- y$crit.consts
#'#The following exists only for step-down direction
#'transformed.pvalues <- y$pval.transf
#'
#'@return
#'For \code{kernel.ADBH.crit}, 
#'a list which elements are:
#'\item{crit.consts}{a vector of critical constants}
#'\item{pval.transf}{a vector of transformed p-values (only for step-down direction).}
#'
#'For \code{kernel.ADBH.fast},
#'a vector of transformed p-values.
#'
#'@name kernel.ADBH
NULL

#'@rdname kernel.ADBH
#'@export
kernel.ADBH.crit <- function(msg = "", stepf, pv.numer, pv.denom, alpha, sorted.pv, bigMem = FALSE, verbose = TRUE){
  l <- length(pv.numer)
  m <- length(stepf)
  # indices of support p-values for which sums are <= alpha * k
  idx <- numeric(m)
  if(bigMem){
    # the small data size allows to use only one matrix
    lst <- tryCatch({
      if(verbose) message(paste(msg, "Using matrix ...", sep = ""))
      # transformed p-values
      y <- numeric(m)
      if(all(pv.numer == pv.denom)){
        # SD case, see (13)
        # fill the matrix
        # columns: indices from 1 to m
        # rows: support p-values
        mat <- sapply(stepf, function(s){ v <- s(pv.numer);return(v / (1 - v)) } )
        # compute the critical constants and transformed p-values
        for(j in 1:l){
          # form the expression in (13)
          temp <- rev(cumsum(sort(mat[j, ], decreasing = TRUE)))
          # critical constants
          idx <- pmax(idx, (temp <= alpha * 1:m) * j)
          # transformed p-values
          k <- which(sorted.pv == pv.numer[j])
          if(length(k)) y[k] <- temp[k]
        }
        # the c(0,pv.numer)[idx+1] is here to allow critical values equal to 0
        lst <- list(crit.consts = c(0,pv.numer)[idx+1], pval.transf = y)
      }else{
        # SU case, see (12)
        # fill the matrix
        # columns: indices from 1 to m
        # rows: support p-values
        mat <- sapply(stepf, function(s){return(s(pv.numer) / (1 - s(pv.denom))) } )
        # compute the critical constants
        for(j in 1:l){
          # form the expression in (12) and deduce critical constants
          idx <- pmax(idx, ( rev(cumsum(sort(mat[j, ], decreasing = TRUE))) <= alpha * 1:m) * j)
        }
        # the c(0,pv.numer)[idx+1] is here to allow critical values equal to 0
        lst <- list(crit.consts = c(0,pv.numer)[idx+1])
      }
      lst
    },
    error = function(err){return(kernel.ADBH.crit(err, stepf, pv.numer, pv.denom, alpha, sorted.pv, FALSE))})
  }
  else{
    # the large data size requires to use chunks
    if(verbose && msg != "") message(paste(msg, "Using matrix with chunking...", sep = ""))
    # size of the chunks
    size <- as.integer(2^25 %/% m)
    # number of chunks
    chunks <- ceiling( l / size)
    # columns: indices from 1 to m
    # rows: support p-values
    mat <- matrix(0, size, m)
    if(all(pv.numer == pv.denom)){
      # SD case, see (13)
      y <- numeric(m)
      for(i in 1:chunks){
        # the min( , l) is here for the last chunk
        chunk <- ((i - 1) * size + 1):min(i * size, l)
        pv.n <- pv.numer[chunk]
        size2 <- length(pv.n)
        # fill the matrix
        for(k in 1:m){
          s_eval <- stepf[[k]](pv.n)
          #size2 may differ from size for the last chunk
          mat[1:size2, k] <- s_eval / (1 - s_eval)
        }
        # compute the critical constants and transformed p-values
        for(j in 1:size2){
          # form the expression in (13)
          temp <- rev(cumsum(sort(mat[j, ], decreasing = TRUE)))
          # critical constants, ((i - 1) * size + j) because we are in chunk number i
          idx <- pmax(idx, (temp <= alpha * 1:m) * ((i - 1) * size + j) )
          # transformed p-values
          k <- which(sorted.pv == pv.n[j])
          if(length(k)) y[k] <- temp[k]
        }
      }
      lst <- list(crit.consts = pv.numer[idx], pval.transf = y)
    }else{
      # SU case, see (12)
      for(i in 1:chunks){
        # the min( , l) is here for the last chunk
        chunk <- ((i - 1) * size + 1):min(i * size, l)
        pv.n <- pv.numer[chunk]
        size2 <- length(pv.n)
        # fill the matrix
        for(k in 1:m){
          #size2 may differ from size for the last chunk
          mat[1:size2, k] <- stepf[[k]](pv.n) / (1 - stepf[[k]](pv.denom))
        }
        # compute the critical constants
        for(j in 1:size2){
          # form the expression in (12) and deduce critical constants, ((i - 1) * size + j) because we are in chunk number i
          idx <- pmax(idx, (rev(cumsum(sort(mat[j, ], decreasing = TRUE))) <= alpha * 1:m) * ((i - 1) * size + j))
        }
      }
      lst <- list(crit.consts = pv.numer[idx])
    }
  }
  return(lst)
}

#'@rdname kernel.ADBH
#'@export
kernel.ADBH.fast <- function(msg = "", stepf, pv.numer, pv.denom, bigMem = FALSE, verbose = TRUE){
  l <- length(pv.numer)
  m <- length(stepf)
  if(bigMem){
    # the small data size allows to use only one matrix
    y <- tryCatch({
      if(verbose) message(paste(msg, "Using matrix ...", sep = ""))
      out <- numeric(l)
      if(all(pv.numer == pv.denom)){
        # SD case, see (13)
        # fill the matrix
        # columns: indices from 1 to m
        # rows: support p-values
        mat <- sapply(stepf, function(s){ v <- s(pv.numer);return(v / (1 - v)) } )
        # remember kernel.ADBH.fast for SD is used with
        # pv.numer equal to the observed p-values,
        # hence l = m
        for(j in 1:l){
          out[j] <- sum(sort(mat[j,], decreasing = TRUE)[1:(m - j + 1)])
        }
      }
      else{
        # SU case, see (12)
        # fill the matrix
        # columns: indices from 1 to m
        # rows: support p-values
        mat <- sapply(stepf, function(s){return(s(pv.numer) / (1 - s(pv.denom))) } )
        # remember kernel.ADBH.fast for SD is used with
        # pv.numer equal to the observed p-values that are <= c.m,
        # hence l <= m
        for(j in 1:l){
          out[j] <- sum(sort(mat[j,], decreasing = TRUE)[1:(m - j + 1)])
        }
      }
      out
    },
    error = function(err){return(kernel.ADBH.fast(err, stepf, pv.numer, pv.denom))})
  }
  else{
    # the large data size requires to use chunks
    if(verbose && msg != "") message(paste(msg, "Using matrix with chunking...", sep = ""))
    # size of the chunks
    size <- as.integer(2^25 %/% m)
    # number of chunks
    chunks <- ceiling(l / size)
    # columns: indices from 1 to m
    # rows: support p-values
    mat <- matrix(0, size, m)
    y <- numeric(l)
    if(all(pv.numer == pv.denom)){
      # SD case, see (13)
      for(i in 1:chunks){
        # the min( , l) is here for the last chunk
        pv.n <- pv.numer[((i - 1) * size + 1):min(i * size, l)]
        size2 <- length(pv.n)
        # fill the matrix
        for(k in 1:m){
          s_eval <- stepf[[k]](pv.n)
          #size2 may differ from size for the last chunk
          mat[1:size2, k] <- s_eval / (1 - s_eval)
        }
        # remember kernel.ADBH.fast for SD is used with
        # pv.numer equal to the observed p-values,
        # hence l = m
        for(j in 1:size2){
          y[(i - 1) * size + j] <- sum((sort(mat[j,], decreasing = TRUE))[1:(m - ((i - 1) * size + j) + 1)])
        }
      }
    }
    else{
      # SU case, see (12)
      for(i in 1:chunks){
        # the min( , l) is here for the last chunk
        pv.n <- pv.numer[((i - 1) * size + 1):min(i * size, l)]
        size2 <- length(pv.n)
        # fill the matrix
        for(k in 1:m){
          #size2 may differ from size for the last chunk
          mat[1:size2, k] <- stepf[[k]](pv.n) / (1 - stepf[[k]](pv.denom))
        }
        # remember kernel.ADBH.fast for SD is used with
        # pv.numer equal to the observed p-values that are <= c.m,
        # hence l <= m
        for(j in 1:size2){
          y[(i - 1) * size + j] <- sum(sort(mat[j,], decreasing = TRUE)[1:(m - ((i - 1) * size + j) + 1)])
        }
      }
    }
  }
  return(y)
}

#'@title [AHSU] and [AHSD] procedures
#'
#'@description
#'Apply the [AHSU] or [AHSD] procedure,
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
#'ADBH.su.fast <- ADBH(pCDFlist, raw.pvalues)
#'ADBH.sd.fast <- ADBH(pCDFlist, raw.pvalues, direction="sd")
#'ADBH.sd.fast$Adjusted
#'ADBH.su.crit <- ADBH(pCDFlist, raw.pvalues, ret.crit.consts=TRUE)
#'ADBH.sd.crit <- ADBH(pCDFlist, raw.pvalues, direction="sd", ret.crit.consts=TRUE)
#'ADBH.sd.crit$Adjusted
#'
#'@templateVar DBR FALSE
#'@template return
#'
#'@name ADBH
#'@export
ADBH <- function(pCDFlist, raw.pvalues, alpha = 0.05, direction = "su", ret.crit.consts = FALSE, bigMem = FALSE, verbose = TRUE){
  m <- length(raw.pvalues)
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
  #        Compute [AHSU] or [AHSD] significant p-values,
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
      # apply the shortcut drawn from Lemma 4, that is
      # c.1 >= the effective critical value associated to min((1-c.m)*alpha/m,c.m)
      pv.list<-short.eff(pv.list.all, min((1-c.m)*alpha/m,c.m) )
      # restrict attention to these values, because c.k needs to be <= c.m
      pv.list <- pv.list[pv.list <= c.m]
      # compute transformed support
      y <- kernel.ADBH.crit("", stepf, pv.list, c.m, alpha, NULL, bigMem, verbose) 
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      # restrict attention to these values, because rejected p-values need to be <= c.m
      obs.pvals<-sorted.pvals[sorted.pvals<=c.m]
      # compute transformed p-values
      y <- kernel.ADBH.fast("", stepf, obs.pvals, c.m, bigMem, verbose) 
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
      # apply the shortcut drawn from Lemma 4, that is
      # c.1 >= the effective critical value associated to alpha/(m + alpha)
      pv.list<-short.eff(pv.list.all, alpha/(m + alpha) )
      # then re-add the observed p-values (needed to compute the adjusted p-values),
      # because we may have removed some of them the shortcut
      pv.list <- sort(unique(c(pv.list, sorted.pvals))) 
      # compute transformed support
      y <- kernel.ADBH.crit("", stepf, pv.list, pv.list, alpha, sorted.pvals, bigMem, verbose) 
      # Find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      # compute transformed p-values
      y <- kernel.ADBH.fast("", stepf, sorted.pvals, sorted.pvals, bigMem, verbose)
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
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, Max.k = m.rej) 
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
  if(ret.crit.consts) output <- c(output, list(Critical.constants = crit.constants)) 
  return(output)
}
