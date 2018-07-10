#'@title Kernel functions for DBR 
#'
#'@description
#'Kernel functions that transform 
#'observed p-values or their support according to [DBR-lambda].
#'The output is used by \code{\link{DBR}}.
#'Additionally, \code{kernel.DBR.crit} computes and returns
#'the critical constants.
#'The end user should not use them.
#'
#'@details
#'When computing critical constants,
#'that is, when using \code{kernel.DBR.crit}, we still need
#'to get transformed p-values to compute
#'the adjusted p-values. Also, note that here the critical constants are computed by
#'the kernel function and not by the principal
#'function \code{\link{DBR}}, contrary to what happens
#'with \code{\link{DBH}}. This is why \code{sorted.pv} is needed.
#'
#'This version: 2018-02-20.
#'
#'@seealso
#'\code{\link{DBR}}, \code{\link{DiscreteFDR}},
#'\code{\link{kernel.DBH}}, \code{\link{DBH}}
#'
#'@templateVar msg TRUE
#'@templateVar stepf TRUE
#'@templateVar pv.numer TRUE
#'@templateVar pv.denom FALSE
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
#'@templateVar lambda TRUE
#'@template param
#'
#'@template example
#'@examples
#'
#'m <- length(raw.pvalues)
#'alpha <- 0.05
#'lambda <- 0.05
#'
#'#Compute the step functions from the supports
#'stepf <- build.stepfuns(pCDFlist)
#'
#'#If not searching for critical constants, we use only the observed p-values
#'sorted.pvals <- sort(raw.pvalues)
#'y <- kernel.DBR.fast("", stepf, sorted.pvals, lambda)
#'
#'#If searching for critical constants, we use (almost) the complete support
#'pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
#'# apply the shortcut drawn from Corollary 3, that is
#'# c.1 >= the effective critical value associated to min((1 - lambda) * alpha/m , lambda)
#'pv.list<-short.eff(pv.list.all, min((1 - lambda) * alpha/m , lambda) )
#'# then re-add the observed p-values (needed to compute the adjusted p-values),
#'# because we may have removed some of them the shortcut
#'pv.list <- sort(unique(c(pv.list, sorted.pvals)))
#'# compute transformed support
#'y <- kernel.DBR.crit("", stepf, pv.list, lambda, alpha, sorted.pvals)
#'crit.constants <- y$crit.consts
#'transformed.pvalues <- y$pval.transf
#'last.index <- y$m.lambda
#'
#'@return
#'For \code{kernel.DBR.crit}, 
#'a list which elements are:
#'\item{crit.consts}{critical constants}
#'\item{pval.transf}{transformed p-values}
#'\item{m.lambda}{last index of observed p-values such that max_i F_i(p) <= lambda,
#'this needs to be passed to \code{\link{DBR}} to compute adjusted p-values).}
#'
#'For \code{kernel.DBR.fast},
#'a vector of transformed p-values.
#'
#'@name kernel.DBR
NULL

#'@rdname kernel.DBR
#'@export
kernel.DBR.crit <- function(msg = "", stepf, pv.numer, lambda, alpha, sorted.pv, bigMem = FALSE, verbose = TRUE){
  l <- length(pv.numer)
  m <- length(stepf)
  # indices of support p-values for which sums are <= alpha * k
  idx <- numeric(m)
  # transformed p-values
  y <- rep(1,m)
  # last index of observed p-values sucht that max_i F_i(p) <= lambda
  k_star <- 0
  if(bigMem){
    # the small data size allows to use only one matrix
    lst <- tryCatch({
      if(verbose) message(paste(msg, "Using matrix ...", sep = ""))
      # fill the matrix
      # columns: indices from 1 to m
      # rows: support p-values
      mat <- matrix(0, l, m)
      mat <- sapply(stepf, function(s){ s(pv.numer) } )
      # screen to keep only values such that max_i F_i(t) <= lambda
      max_F <- apply(mat, 1, max)
      idx_lambda <- which(max_F <= lambda)
      for(j in idx_lambda){
        # work only on values such that max_i F_i(t) <= lambda
        temp <- rev(cumsum(sort(mat[j,], decreasing = TRUE))) / (1 - lambda)
        # critical constants
        idx <- pmax(idx, (temp <= alpha * 1:m) * j)
        # transformed p-values
        k <- which(sorted.pv == pv.numer[j])
        if(length(k)) {
          k_star<-k
          y[k] <- temp[k]
        }
      }
      # the c(0,pv.numer)[idx+1] is here to allow critical values equal to 0
      lst <- list(crit.consts = c(0,pv.numer)[idx+1], pval.transf = y, m.lambda = k_star)
    },
    error = function(err){return(kernel.DBR.crit(err, stepf, pv.numer, lambda, alpha, sorted.pv))})
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
    for(i in 1:chunks){
      # the min( , l) is here for the last chunk
      pv.n <- pv.numer[((i - 1) * size + 1):min(i * size, l)]
      size2 <- length(pv.n)
      # fill the matrix
      for(j in 1:m){
        #size2 may differ from size for the last chunk
        mat[1:size2, j] <- stepf[[j]](pv.n)
      }
      # compute the critical constants and transformed p-values
      # screen to keep only values such that max_i F_i(t) <= lambda
      max_F <- apply(mat[1:size2,], 1, max)
      idx_lambda <- which(max_F <= lambda)
      for(j in idx_lambda){
        temp <- rev(cumsum(sort(mat[j, ], decreasing = TRUE))) / (1 - lambda)
        # critical constants, ((i - 1) * size + j) because we are in chunk number i
        idx <- pmax(idx, (temp <= alpha * 1:m) * ((i - 1) * size + j))
        # transformed p-values
        k <- which(sorted.pv == pv.n[j])
        if(length(k)) {
          k_star<-k
          y[k] <- temp[k]
        }
      }
      if(length(idx_lambda) < size2){
        # we have screened through all values such that max_i F_i(t) <= lambda,
        # so we break the loop over chunks to prevent useless computations
        break
      }
    }
    # the c(0,pv.numer)[idx+1] is here to allow critical values equal to 0
    lst <- list(crit.consts = c(0,pv.numer)[idx+1], pval.transf = y, m.lambda = k_star )
  }
  return(lst)
}

#'@rdname kernel.DBR
#'@export
kernel.DBR.fast <- function(msg = "", stepf, pv.numer, lambda, bigMem = FALSE, verbose = TRUE){
  m <- length(stepf)
  if(bigMem){
    # the small data size allows to use only one matrix
    y <- tryCatch({
      if(verbose) message(paste(msg, "Using matrix ...", sep = ""))
      out <- rep(1, m)
      # fill the matrix
      # columns: indices from 1 to m
      # rows: support p-values
      mat <- sapply(stepf, function(s){ s(pv.numer) } )
      # screen to keep only p-values such that max_i F_i(p) <= lambda
      max_F <- apply(mat, 1, max)
      idx_lambda <- which(max_F <= lambda)
      for(j in idx_lambda){
        out[j] <- sum(sort(mat[j,], decreasing = TRUE)[1:(m - j + 1)]) / ((1 - lambda) * j)
      }
      out
    },
    error = function(err){return(kernel.DBR.fast(err, stepf, pv.numer, lambda))})
  }
  else{
    # the large data size requires to use chunks
    if(verbose && msg != "") message(paste(msg, "Using matrix with chunking...", sep = ""))
    y <- rep(1, m)
    # size of the chunks
    size <- as.integer(2^25 %/% m)
    # number of chunks
    chunks <- ceiling(m / size)
    # columns: indices from 1 to m
    # rows: support p-values
    mat <- matrix(0, size, m)
    for(i in 1L:chunks){
      # the min( , m) is here for the last chunk
      pv.n <- pv.numer[((i - 1) * size + 1):min(i * size, m)]
      size2 <- length(pv.n)
      # fill the matrix
      for(k in 1:m){
        #size2 may differ from size for the last chunk
        mat[1:size2, k] <- stepf[[k]](pv.n)
      }
      # screen to keep only values such that max_i F_i(t) <= lambda
      max_F <- apply(mat[1:size2,], 1, max)
      idx_lambda <- which(max_F <= lambda)
      m_lambda <- length(idx_lambda)
      for(j in idx_lambda){
        # ((i - 1) * size + j) because we are in chunk number i
        y[(i - 1) * size + j] <- sum(sort(mat[j,], decreasing = TRUE)[1:(m - ((i - 1) * size + j) + 1)]) / ((1 - lambda) * ((i - 1) * size + j))
      }
      if(length(idx_lambda) < size2){
        # we have screened through all values such that max_i F_i(t) <= lambda,
        # so we break the loop over chunks to prevent useless computations
        break
      }
    }
  }
  return(y)
}

#'@title [DBR-lambda] procedure
#'
#'@description
#'Apply the [DBR-lambda] procedure,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'This version: 2018-02-19.
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
#'DBR.su.fast <- DBR(pCDFlist, raw.pvalues)
#'DBR.su.crit <- DBR(pCDFlist, raw.pvalues, ret.crit.consts=TRUE)
#'
#'@templateVar DBR TRUE
#'@template return
#'
#'@name DBR
#'@export
DBR <- function(pCDFlist, raw.pvalues, alpha = 0.05, lambda = 0.05, ret.crit.consts = FALSE, bigMem = FALSE, verbose = FALSE){
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
    pv.list<-short.eff(pv.list.all, min((1 - lambda) * alpha/m , lambda) )
    # then re-add the observed p-values (needed to compute the adjusted p-values),
    # because we may have removed some of them the shortcut
    pv.list <- sort(unique(c(pv.list, sorted.pvals)))
    # compute transformed support
    y <- kernel.DBR.crit("", stepf, pv.list, lambda, alpha, sorted.pvals, bigMem, verbose)
    # find critical constants
    crit.constants<-y$crit.consts
    idx <- which(sorted.pvals <= crit.constants)   
  }
  else{
    # compute transformed p-values
    y <- kernel.DBR.fast("", stepf, sorted.pvals, lambda, bigMem, verbose)
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
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, Max.k = m.rej, Lambda = lambda)
  if(ret.crit.consts){
    output <- c(output, list(Critical.constants = crit.constants))
    # compute adjusted p-values
    # recall that transformed p-values where max_i F_i(p) > lambda,
    # that is for indices > y$m.lambda, are set to 1
    pv.adj <- rev(cummin(rev(  pmin(y$pval.transf / c(seq_len(y$m.lambda), rep(1,m-y$m.lambda) ) , 1)  )))
  }
  else{
    # compute adjusted p-values
    pv.adj <- rev(cummin(rev( pmin(y,1) )))
  }
  # add adjusted p-values to output list
  ro <- order(o)
  output <- c(output, list(Adjusted = pv.adj[ro]))
  return(output)
}