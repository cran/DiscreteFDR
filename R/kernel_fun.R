#'@title Kernel functions for discrete.BH
#'
#'@description
#'Kernel functions that transform 
#'observed p-values or their support according to [HSU], [HSD],
#'[AHSU] and [AHSD]. 
#'The output is used by \code{\link{discrete.BH}}.
#'Additionally, \code{kernel.DBH.crit} and \code{kernel.ADBH.crit}
#'compute and return the critical constants.
#'The end user should not use them.
#'
#'@details
#'When computing critical constants under step-down, that is,
#'when using \code{kernel.DBH.crit} or \code{kernel.ADBH.crit}
#'with \code{pv.numer = pv.denom} (i.e. the step-down case),
#'we still need to get transformed p-values to compute
#'the adjusted p-values.
#'
#'This version: 2018-11-13.
#'
#'@seealso
#'\code{\link{discrete.BH}}, \code{\link{DiscreteFDR}},
#'\code{\link{DBR}}, \code{\link{kernel.DBR}}
#'
#'@templateVar stepf TRUE
#'@templateVar pv.numer TRUE
#'@templateVar pv.denom TRUE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv TRUE
#'@templateVar pCDFlist FALSE
#'@templateVar raw.pvalues FALSE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@templateVar adaptive FALSE
#'@template param 
#'
#'@template example
#'@examples
#'
#'m <- length(raw.pvalues)
#'alpha <- 0.05
#'
#'# Compute the step functions from the supports
#'stepf <- build.stepfuns(pCDFlist)
#'
#'# We stay in a step-down context, where pv.numer = pv.denom,
#'# for the sake of simplicity
#'
#'# If not searching for critical constants, we use only the observed p-values
#'sorted.pvals <- sort(raw.pvalues)
#'y.DBH.fast <- kernel.DBH.fast(stepf, sorted.pvals, sorted.pvals)
#'y.ADBH.fast <- kernel.ADBH.fast(stepf, sorted.pvals, sorted.pvals)
#'
#'# If searching for critical constants, we use (almost) the complete support
#'pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
#'# apply the shortcut drawn from Lemma 3, that is
#'# c.1 >= the effective critical value associated to (alpha/m)/(1 + alpha/m)
#'pv.list <- short.eff(pv.list.all, alpha/(m + alpha))
#'# then re-add the observed p-values (needed to compute the adjusted p-values),
#'# because we may have removed some of them by the shortcut
#'pv.list <- sort(unique(c(pv.list, sorted.pvals)))
#'# compute transformed support
#'y.DBH.crit <- kernel.DBH.crit(stepf, pv.list, pv.list, alpha, sorted.pvals)
#'y.ADBH.crit <- kernel.ADBH.crit(stepf, pv.list, pv.list, alpha, sorted.pvals)
#'# critical constants
#'y.DBH.crit$crit.consts
#'y.ADBH.crit$crit.consts
#'# The following exist only for step-down direction
#'y.DBH.crit$pval.transf
#'y.ADBH.crit$pval.transf
#'
#'@return
#'For \code{kernel.DBH.fast} and \code{kernel.ADBH.fast}, a vector of transformed p-values is returned.
#'\code{kernel.DBH.crit} and \code{kernel.ADBH.crit} return a list object with critical constants (\code{$crit.consts})
#'and transformed p-values (\code{$pval.transf}), but only if \code{pv.numer = pv.denom}.
#'
#'@name kernel
NULL

#'@rdname kernel
#'@export
kernel.DBH.fast <- function(stepf, pv.numer, pv.denom){
  m <- length(stepf)
  # the output is the vector y= \sum_{i=1}^m F_i(pv.numer)/(1 - F_i(pv.denom))
  # note that pv.number is either all values of all supports of the p-values
  # (when searching for the critical values), or only the observed p-values
  # (when skipping the critical values computation)
  
  # possibly large data size requires to use a for loop
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
  return(y)
}

#'@rdname kernel
#'@export
kernel.DBH.crit <- function(stepf, pv.numer, pv.denom, alpha, sorted.pv){
  m <- length(stepf)
  y <- kernel.DBH.fast(stepf, pv.numer, pv.denom)
  crit <- sapply(1:m, function(i){pv.numer[which.max(y[y <= i * alpha])]})
  lst <- list(crit.consts = crit)
  if(all(pv.numer == pv.denom)){
    # SD case
    # get the indices of the observed p-values inside pv.list
    idx <- sapply(1:m, function(i){which(pv.numer == sorted.pv[i])})
    lst$pval.transf = y[idx]
  }
  
  return(lst)
}

#'@rdname kernel
#'@export
kernel.ADBH.fast <- function(stepf, pv.numer, pv.denom){
  l <- length(pv.numer)
  m <- length(stepf)
  
  # possibly large data size requires to use chunks
  # size of the chunks
  size <- as.integer(2^25 %/% m) # 2^25 numeric elements require ~256 MiB of RAM; best size for most performance
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
  return(y)
}

#'@rdname kernel
#'@export
kernel.ADBH.crit <- function(stepf, pv.numer, pv.denom, alpha, sorted.pv){
  l <- length(pv.numer)
  m <- length(stepf)
  # indices of support p-values for which sums are <= alpha * k
  idx <- numeric(m)

  # possibly large data size requires to use chunks
  # size of the chunks
  size <- as.integer(2^25 %/% m) # 2^25 numeric elements require ~256 MiB of RAM; best size for most performance
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
  return(lst)
}

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
#'When computing critical constants, that is,
#'when using \code{kernel.DBR.crit}, we still need
#'to get transformed p-values to compute the adjusted 
#'p-values.
#'
#'This version: 2018-11-13.
#'
#'@seealso
#'\code{\link{DBR}}, \code{\link{DiscreteFDR}},
#'\code{\link{kernel}}, \code{\link{discrete.BH}}
#'
#'@templateVar stepf TRUE
#'@templateVar pv.numer TRUE
#'@templateVar pv.denom FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv TRUE
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
#'y <- kernel.DBR.fast(stepf, sorted.pvals, lambda)
#'
#'#If searching for critical constants, we use (almost) the complete support
#'pv.list.all <- unique(sort(as.numeric(unlist(pCDFlist))))
#'# apply the shortcut drawn from Corollary 3, that is
#'# c.1 >= the effective critical value associated to min((1 - lambda) * alpha/m , lambda)
#'pv.list <- short.eff(pv.list.all, min((1 - lambda) * alpha/m , lambda))
#'# then re-add the observed p-values (needed to compute the adjusted p-values),
#'# because we may have removed some of them the shortcut
#'pv.list <- sort(unique(c(pv.list, sorted.pvals)))
#'# compute transformed support
#'y <- kernel.DBR.crit(stepf, pv.list, lambda, alpha, sorted.pvals)
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
kernel.DBR.fast <- function(stepf, pv.numer, lambda){
  m <- length(stepf)
  # possibly large data size requires to use chunks
  y <- rep(1, m)
  # size of the chunks
  size <- as.integer(2^25 %/% m) # 2^25 numeric elements require ~256 MiB of RAM; best size for most performance
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
  return(y)
}

#'@rdname kernel.DBR
#'@export
kernel.DBR.crit <- function(stepf, pv.numer, lambda, alpha, sorted.pv){
  l <- length(pv.numer)
  m <- length(stepf)
  # indices of support p-values for which sums are <= alpha * k
  idx <- numeric(m)
  # transformed p-values
  y <- rep(1,m)
  # last index of observed p-values sucht that max_i F_i(p) <= lambda
  k_star <- 0
  
  # possibly large data size requires to use chunks
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
  lst <- list(crit.consts = c(0, pv.numer)[idx + 1], pval.transf = y, m.lambda = k_star)
  return(lst)
}