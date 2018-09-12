#'@description
#' This package implements the [HSU], [HSD],
#' [AHSU], [AHSD] and [HBR-lambda] procedures for
#' discrete tests (see References). 
#'
#'@details
#' The functions are reorganised from the reference paper in the following way.
#' \code{\link{DBH}} (for Discrete Benjamini-Hochberg) implements
#' [HSU] and [HSD], \code{\link{ADBH}} (the "A" stands for Adaptive) implements [AHSU] and [AHSD], 
#' and \code{\link{DBR}} (for Discrete Blanchard-Roquain) implements [HBR-lambda].
#' Their main arguments are a vector
#' of raw observed p-values, and a list
#' of the same length, which elements are the discrete supports
#' of the CDFs of the p-values.
#' 
#' The function \code{\link{fisher.pvalues.support}} allows to compute
#' such p-values and support in the framework of a Fisher's
#' exact test of association. It has been inspired by an help
#' page of the package \code{\link[discreteMTP]{discreteMTP}}.
#' 
#' The function \code{\link{fast.Discrete}} is a wrapper for \code{\link{fisher.pvalues.support}}
#' and  \code{\link{DBH}} (or  \code{\link{ADBH}}) which allows to apply discrete procedures
#' directly to a data set of contingency tables.
#' 
#' We also provide the \code{\link{amnesia}} data set, used in 
#' our examples and in our paper. It is basically the \code{\link[discreteMTP]{amnesia}} data set
#' of package \code{\link[discreteMTP]{discreteMTP}}, but slightly reformatted (the difference lies in column 3).
#' 
#' No other function of the package should
#' be used, they are only internal functions called
#' by the main ones.
#' 
#' @section References:
#' S. Döhler, G. Durand and E. Roquain (2018). New FDR bounds for discrete and heterogeneous tests. Electronic Journal of Statistics, Volume 12, Number 1 (2018) \href{https://projecteuclid.org/euclid.ejs/1528855551}{link}.
"_PACKAGE"