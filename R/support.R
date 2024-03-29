#' @title
#' Computing discrete p-values and their support for binomial and Fisher's
#' exact tests
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Computes discrete raw p-values and their support
#' for binomial test or Fisher's exact test applied to 2x2 contingency tables
#' summarizing counts coming from two categorical measurements.
#'
#' **Note**: In future versions, this function will be removed. Generation of
#' p-value supports for different exact tests, including Fisher's, will be
#' moved to a separate package.
#' 
#' @details
#' Assume that each contingency tables compares two variables and resumes the
#' counts of association or not with a condition. This can be resumed in the
#' following table:
#' \tabular{lccc}{
#' \tab Association \tab No association  \tab      Total      \cr
#'      Variable 1  \tab    \eqn{X_1}    \tab    \eqn{Y_1}    \tab \eqn{N_1} \cr
#'      Variable 2  \tab    \eqn{X_2}    \tab    \eqn{Y_2}    \tab \eqn{N_2} \cr
#'      Total       \tab \eqn{X_1 + X_2} \tab \eqn{Y_1 + Y_2} \tab \eqn{N_1 + N_2}
#' }
#' If `input="noassoc"`, `counts` has four columns which respectively contain,
#' \eqn{X_1}, \eqn{Y_1}, \eqn{X_2} and \eqn{Y_2}. If `input="marginal"`,
#' `counts` has four columns which respectively contain \eqn{X_1}, \eqn{N_1},
#' \eqn{X_2} and \eqn{N_2}.
#' 
#' If `input="HG2011"`, we are in the situation of the [amnesia] data set as
#' in Heller & Gur (2011, see References). Each contingency table is obtained
#' from one variable which is compared to all other variables of the study. That
#' is, counts for "second variable" are replaced by the sum of the counts of the
#' other variables:
#' \tabular{lccc}{
#' \tab Association            \tab No association            \tab Total                     \cr
#'      Variable \eqn{j}       \tab \eqn{X_j}                 \tab \eqn{Y_j}                 \tab \eqn{N_j} \cr
#'      Variables \eqn{\neq j} \tab \eqn{\sum_{i \neq j} X_i} \tab \eqn{\sum_{i \neq j} Y_i} \tab \eqn{\sum_{i \neq j} N_i} \cr
#'      Total                  \tab \eqn{\sum X_i}            \tab \eqn{\sum Y_i}            \tab \eqn{\sum N_i}
#' }
#' Hence `counts` needs to have only two columns which respectively contain \eqn{X_j} and \eqn{Y_j}.
#'
#' The code for the computation of the p-values of Fisher's exact test is
#' inspired by the example in the help page of `p.discrete.adjust` of package
#' `discreteMTP`, which is no longer available on CRAN.
#'
#' See the Wikipedia article about Fisher's exact test, paragraph Example, for
#' a good depiction of what the code does for each possible value of
#' `alternative`.
#'
#' @seealso
#' [fisher.test]
#' 
#' @param counts        a data frame of two or four columns and any number of
#'                      lines; each line represents a 2x2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument, see
#'                      Details.
#' @param alternative   same argument as in [stats::fisher.test]. The three
#'                      possible values are `"greater"` (default),
#'                      `"two.sided"` or `"less"` and you can specify
#'                      just the initial letter.
#' @param input         the format of the input data frame, see Details. The
#'                      three possible values are `"noassoc"` (default),
#'                      `"marginal"` or `"HG2011"` and you can specify
#'                      just the initial letter.
#' 
#' @template example
#' @template exampleHG
#' 
#' @return
#' A list of two elements:
#' \item{raw}{raw discrete p-values.}
#' \item{support}{a list of the supports of the CDFs of the p-values.
#' Each support is represented by a vector in increasing order.}
#' 
#' @references
#' R. Heller and H. Gur (2011). False discovery rate controlling procedures for
#'   discrete tests. arXiv preprint.
#'   [arXiv:1112.4627v2](https://arxiv.org/abs/1112.4627v2).
#'
#' "Fisher's exact test", Wikipedia, The Free Encyclopedia, accessed 2018-03-20,
#' [link](https://en.wikipedia.org/w/index.php?title=Fisher's_exact_test&oldid=823327889).
#' 
#' @importFrom stats dhyper phyper pbinom
#' @importFrom lifecycle deprecate_soft
#' @export
fisher.pvalues.support <- function(counts, alternative = "greater", input = "noassoc"){
  deprecate_soft("1.3.7", "fast.Discrete()",
                 details = paste("In future versions of this package, this",
                                 "function will be removed. Generation of",
                                 "p-value supports for different exact tests,",
                                 "including Fisher's, will be moved to a",
                                 "separate package."))
  
  input <- match.arg(input, c("noassoc", "marginal", "HG2011"))
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  
  number.items <- nrow(counts)
  switch (input, 
          noassoc = {
            A11 <- counts[, 1]
            A12 <- counts[, 2]
            n <- A11 + A12
            #A21 <- counts[, 3]
            #A22 <- counts[, 4]
            A1. <- A11 + counts[, 3]
            A2. <- A12 + counts[, 4]
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the j-th test of no association
            # between two variables and a condition :  
            #                   Association    No association   Marginal counts
            #  Variable 1       A11[j]         A12[j]           n[j]
            #  Variable 2       A21[j]         A22[j]           N[j]-n[j]
            #  All variables    A1.[j]         A2.[j]           N[j]
          },
          marginal = {
            A11 <- counts[, 1]
            n <- counts[, 2]
            A12 <- n - A11
            #A21 <- counts[, 3]
            #A22 <- counts[, 4] - A21
            A1. <- A11 + counts[, 3]
            A2. <- counts[, 4] + n - A1.
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the j-th test of no association
            # between two variables and a condition :  
            #                   Association    No association   Marginal counts
            #  Variable 1       A11[j]         A12[j]           n[j]
            #  Variable 2       A21[j]         A22[j]           N[j]-n[j]
            #  All variables    A1.[j]         A2.[j]           N[j]
          },
          HG2011 = {
            A11 <- counts[, 1]
            #A21 <- sum(counts[, 1]) - A11
            A12 <- counts[, 2]
            #A22 <- sum(counts[, 2]) - A12
            A1. <- rep(sum(counts[, 1]), number.items)
            A2. <- rep(sum(counts[, 2]), number.items)
            n <- A11 + A12
            #N <- A1. + A2.
            # Entry j in each of the above vectors is a count for the test of no association
            # between variable j and the condition, compared to all other variables :  
            #                    Association    No association   Marginal counts
            #  Variable j        A11[j]         A12[j]           n[j]
            #  Other variables   A21[j]         A22[j]           N[j]-n[j]
            #  All variables     A1.[j]         A2.[j]           N[j]
          }
  )
  
  pCDFlist <- vector("list", number.items)
  raw.pvalues <- numeric(number.items)
  
  k <- pmin(n, A1.)
  l <- pmax(0, n - A2.)
  switch(alternative, greater =
           for (i in 1:number.items){
             x <- l[i]:k[i]
             # the "-1" below is because lower.tail = FALSE computes P[X > x],
             # and we want P[X >= x]=P[X > x-1]
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, phyper(x-1, A1.[i], A2.[i], n[i], lower.tail = FALSE)))
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i] - l[i] + 1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- rev(pCDFlist[[i]])
           },
         less =
           for (i in 1:number.items){
             x <- l[i]:k[i]
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, phyper(x, A1.[i], A2.[i], n[i], lower.tail = TRUE)))
             # the "+1" below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i] - l[i] + 1]
           },
         two.sided =
           for (i in 1:number.items){
             x <- l[i]:k[i]
             atoms <- dhyper(x, A1.[i], A2.[i], n[i])
             # ensure that probabilities sum up to 1 (sometimes, needs to be done multiple times)
             newsum <- sum(atoms)
             while(newsum < 1){
               oldsum <- newsum
               atoms <- atoms/newsum
               newsum <- sum(atoms)
               if(oldsum == newsum) break;
             }
			       # pmin/pmax below is to account for machine rounding issues
             pCDFlist[[i]] <- pmax(0, pmin(1, sapply(x, function(nu){sum(atoms[which(atoms <= atoms[nu + 1 - l[i]])])})))
             # the "+1" above and below is because vectors start with index 1 in R: x[A11[i]+1]=A11[i]
             raw.pvalues[i] <- pCDFlist[[i]][A11[i] - l[i] + 1]
             # we want to have pCDFlist[[i]] in increasing order:
             pCDFlist[[i]] <- sort(pCDFlist[[i]])
           }
  )
  
  return(list(raw = raw.pvalues, support = pCDFlist))
}
