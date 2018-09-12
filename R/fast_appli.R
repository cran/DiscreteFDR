#'@title Fast application of discrete procedures
#'
#'@description
#'Apply the [HSU], [HSD], [AHSU] or [AHSD] procedure,
#'without computing the critical constants,
#'to a data set of 2 x 2 contingency tables using Fisher's exact tests.
#'
#'@details
#'This version: 2018-08-21.
#'
#'@param counts a data frame of 2 or 4 columns and any number of lines,
#'each line representing a 2 x 2 contingency table to test.
#'The number of columns and what they must contain depend on the value of the \code{input} argument, 
#'see Details of \code{\link{fisher.pvalues.support}}.
#'@param alternative same argument as in \code{\link{fisher.test}}. The three possible values are \code{"greater"} 
#'(default), \code{"two.sided"} or \code{"less"} and you can specify just the initial letter.
#'@param input the format of the input data frame, see Details of \code{\link{fisher.pvalues.support}}.
#'The three possible values are \code{"noassoc"} 
#'(default), \code{"marginal"} or \code{"HG2011"} and you can specify just the initial letter.
#'@templateVar msg FALSE
#'@templateVar stepf FALSE
#'@templateVar pv.numer FALSE
#'@templateVar pv.denom FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted.pv FALSE
#'@templateVar bigMem TRUE
#'@templateVar verbose TRUE
#'@templateVar pCDFlist FALSE
#'@templateVar raw.pvalues FALSE
#'@templateVar direction TRUE
#'@templateVar ret.crit.consts FALSE
#'@templateVar sorted.num FALSE
#'@templateVar t FALSE
#'@templateVar lambda FALSE
#'@template param
#'@param adaptive a boolean specifying whether to conduct an adaptive procedure or not.
#'
#'@examples
#'
#'X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#'X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#'N1 <- rep(148, 9)
#'N2 <- rep(132, 9)
#'Y1 <- N1-X1
#'Y2 <- N2-X2
#'df <- data.frame(X1,Y1,X2,Y2)
#'df
#'
#'DBH.su <- fast.Discrete(counts=df, input="noassoc", direction="su")
#'DBH.sd <- fast.Discrete(counts=df, input="noassoc", direction="sd")
#'DBH.sd$Adjusted
#'ADBH.su <- fast.Discrete(counts=df, input="noassoc", direction="su", adaptive=TRUE)
#'ADBH.sd <- fast.Discrete(counts=df, input="noassoc", direction="sd", adaptive=TRUE)
#'ADBH.sd$Adjusted
#'
#'@return 
#'A list whose elements are:
#'\item{Rejected}{rejected raw p-values}
#'\item{Indices}{indices of rejected hypotheses}
#'\item{Max.k}{number of rejections}
#'\item{Alpha}{maximum significance level for which a rejection occured, that is \eqn{Alpha = alpha * Max.k / m}}
#'\item{Adjusted}{adjusted p-values (only for step-down direction).}
#'
#'@name fast.Discrete
#'@export
fast.Discrete <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, direction = "su", adaptive = FALSE, bigMem = FALSE, verbose = FALSE){
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  if (adaptive){
    return(ADBH(pCDFlist, raw.pvalues, alpha, direction, FALSE, bigMem, verbose))
  }else{
    return(DBH(pCDFlist, raw.pvalues, alpha, direction, FALSE, bigMem, verbose))
  }
}