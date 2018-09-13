# Purpose: Package documentation
# Updated: 180912

#' @useDynLib RNOmni
#' @importFrom Rcpp sourceCpp
NULL

#' RNOmni: Rank-Normal Omnibus Association Test
#' 
#' Implementation of genetic association tests that utilize the rank-based inverse
#' normal transformation (INT). The main function is \code{\link{RNOmni}}, which performs
#' an omnibus, INT-based association test. The omnibus test combines the direct \code{\link{DINT}}
#' and indirect \code{\link{IINT}} INT-based tests. The transformation itself is described 
#' under \code{\link{rankNorm}}. See the vignette for additional details. 
#' 
#' @author Zachary R. McCaw
#' @docType package
#' @name RNOmni-help
NULL