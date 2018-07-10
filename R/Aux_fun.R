#'@title Matching raw p-values with supports 
#'
#'@description
#'Constructs the observed p-values from the raw observed p-values, 
#'by rounding them to their nearest neighbour matching
#'with the supports of their respective
#'CDFs (as in function \code{\link[discreteMTP]{p.discrete.adjust}} 
#'of package \code{\link[discreteMTP]{discreteMTP}}).
#'The end user should not use it.
#'
#'@details
#'Well computed raw p-values should already 
#'belong to their respective CDF support.
#'So this function is called at the beginning
#'of \code{\link{DBH}}, \code{\link{ADBH}},
#'and \code{\link{DBR}}, just in case
#'raw p-values are biased.
#'
#'For each raw p-value that needs to be rounded,
#'a warning is issued.
#'
#'This version: 2017-08-16.
#'
#'@seealso
#'\code{\link{DBH}}, \code{\link{ADBH}}, \code{\link{DBR}}
#'
#'@templateVar msg FALSE
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha FALSE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem FALSE
#'@templateVar verbose FALSE
#'@templateVar pCDFlist TRUE
#'@templateVar raw.pvalues TRUE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@template param
#'
#'@examples
#'toyList<-list(c(0.3,0.7,1),c(0.1,0.65,1))
#'toyRaw1<-c(0.3,0.65)
#'match.pvals(toyList,toyRaw1)
#'toyRaw2<-c(0.31,0.6)
#'match.pvals(toyList,toyRaw2)
#'
#'@return
#'A vector where each raw p-value has been
#'replaced by its nearest neighbour.
#'
#'@name match.pvals
#'@export
match.pvals <- function(pCDFlist, raw.pvalues){
  m <- length(raw.pvalues)
  pvec <- raw.pvalues
  in.CDF <- numeric(m)
  for (k in (1:m)) {
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
}

#'@title Building step functions from \code{pCDFlist} 
#'
#'@description
#'Creates a list of step functions out of p-value
#'CDF supports. That is, creates the list of the CDFs.
#'The end user should not use it.
#'
#'@details
#'The functions returned are the ones denoted by \eqn{F_i} in the
#'reference paper, see \code{\link{DiscreteFDR}}.
#'
#'This version: 2017-09-09.
#'
#'@templateVar msg FALSE
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha FALSE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem FALSE
#'@templateVar verbose FALSE
#'@templateVar pCDFlist TRUE
#'@templateVar raw.pvalues FALSE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@template param
#'
#'@examples
#'toyList<-list(c(0.3,0.7,1),c(0.1,0.65,1))
#'toyStep<-build.stepfuns(toyList)
#'toyStep[[1]](0.2)
#'toyStep[[2]](0.2)
#'toyStep[[1]](0.65)
#'toyStep[[2]](0.65)
#'
#'@return
#'A list of CDFs.
#'
#'@importFrom stats stepfun
#'@name build.stepfuns
#'@export
build.stepfuns <- function(pCDFlist){
  m <- length(pCDFlist)
  stepf <- vector("list", m)
  for(l in 1:m){
    stepf[[l]] <- stepfun(pCDFlist[[l]] , c(0, pCDFlist[[l]]))
  }
  return(stepf)
}

#'@title Shortcuts for critical values computation
#'
#'@description
#'Extracts all values from a sorted vector
#'that are greater than or equal to the effective
#'critical value associated to a threshold.
#'
#'@details
#'The effective critical value associated to \code{t} is
#'the largest value of \code{sorted.num} that is less than 
#'or equal to \code{t}.
#'
#'This version: 2018-02-12.
#'
#'@templateVar msg FALSE
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha FALSE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem FALSE
#'@templateVar verbose FALSE
#'@templateVar pCDFlist FALSE
#'@templateVar raw.pvalues FALSE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num TRUE
#'@templateVar t TRUE
#'@templateVar lambda FALSE
#'@template param
#'
#'@examples
#'x<-c(0.1,0.2,0.3,0.4)
#'short.eff(x,0.2)
#'short.eff(x,0.35)
#'
#'@return
#'A subvector of \code{sorted.num}.
#'
#'@name short.eff
#'@export
short.eff <- function(sorted.num,t){
  w<-which(sorted.num > t )
  # sorted_num is sorted in ascending order,
  # so the effective critical value associated
  # to t is indexed by w[1]-1, and the subvector
  # we are searching for is indexed by c(w[1]-1,w)
  return(sorted.num[c(w[1]-1,w)])
}

#'@title Computing discrete p-values and their support for Fisher's exact tests
#'
#'@description
#'Computes discrete raw p-values and their support
#'for the test of no association between two categorical 
#'variables in 2 x 2 contingency tables using Fisher's exact tests.
#'
#'@details
#'The code for this function is inspired 
#'from the example in the help page
#'of \code{\link[discreteMTP]{p.discrete.adjust}}.
#'
#'See the Wikipedia article about Fisher's
#'exact test, paragraph Example, for a good depiction
#'of what the code does for each possible value
#'of \code{alternative}.
#'
#'This version: 2018-03-20.
#'
#'@seealso
#'\code{\link[discreteMTP]{p.discrete.adjust}}, \code{\link{fisher.test}}
#'
#'@param counts a data frame of 3 columns and any number of lines,
#'each line being an item for which we want to perform a test.
#'The first column is the name of the item,
#'the second is the count of associations between the item and the condition,
#'the third is the count of no associations.
#'@param alternative same argument as in \code{\link{fisher.test}}.
#'
#'@template example
#'
#'@return
#'A list of two elements:
#'\item{raw}{raw discrete p-values}
#'\item{support}{a list of the supports of the CDFs of the p-values.
#'Each support is represented by a vector in increasing order.}
#'
#'@section References:
#' "Fisher's exact test", Wikipedia, The Free Encyclopedia,
#' accessed 2018-03-20,
#' \href{https://en.wikipedia.org/w/index.php?title=Fisher\%27s_exact_test&oldid=823327889}{link}
#'
#'@importFrom stats dhyper phyper
#'@name fisher.pvalues.support
#'@export
fisher.pvalues.support<-function(counts, alternative = "greater"){
  A11 <- counts[,2]
  #A21 <- sum(counts[,2]) - A11
  A12 <- counts[,3]
  #A22 <- sum(counts[,3]) - A12
  A1. <- sum(counts[,2])
  A2. <- sum(counts[,3])   
  n <- A11 + A12
  #N<- A1. + A2.
  # Entry j in each of the four vectors is the data for the test of no association
  # between item j and the condition :  
  #                        Item j    Other items  All items
  #  Condition             A11[j]    A21[j]       A1.     
  #  No condition          A12[j]    A22[j]       A2.
  #  Total                 n[j]      N-n[j]       N		
  k <- pmin(n,A1.)
  pCDFlist <- list()
  number.items <- nrow(counts)
  raw.pvalues <- numeric(number.items)
  alternative <- char.expand(alternative, c("two.sided", "less", "greater"))
  switch(alternative,greater=
           for (i in 1:number.items){
             x <- 0:k[i]
             # the "-1" below is because lower.tail = FALSE computes P[X > x],
             # and we want P[X >= x]=P[X > x-1]
             pCDFlist[[i]] <- phyper(x-1, A1., A2., n[i], lower.tail = FALSE)
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i]+1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- rev(pCDFlist[[i]])
           },
         less=
           for (i in 1:number.items){
             x <- 0:k[i]
             pCDFlist[[i]] <- phyper(x, A1., A2., n[i], lower.tail = TRUE)
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i]+1]
           },
         two.sided=
           for (i in 1:number.items){
             x <- 0:k[i]
             atoms <- dhyper(x, A1., A2., n[i])
             pCDFlist[[i]] <- sapply(x, function(nu){sum(atoms[which(atoms<=atoms[nu+1])])})
             # the "+1" above and below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i]+1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- sort(pCDFlist[[i]])
           }
  )
  return(list(raw = raw.pvalues, support = pCDFlist))
}
