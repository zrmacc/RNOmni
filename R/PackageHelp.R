# Purpose: Package documentation
# Updated: 19/01/11

#' @useDynLib RNOmni
#' @importFrom Rcpp sourceCpp
NULL

#' RNOmni: Rank-Normal Omnibus Association Testing
#' 
#' Implementation of genetic association tests that incorporate the rank-based inverse
#' normal transformation (INT) \code{\link{rankNorm}}. The direct-INT \code{\link{DINT}} test directly
#' transforms the outcome, whereas the indirect-INT \code{\link{IINT}} test forms residuals prior to transformation. 
#' The omnibus INT \code{\link{OINT}} test adaptively combines the D-INT and I-INT tests into a single
#' robust and statistically powerful procedure.
#' 
#' @author Zachary R. McCaw
#' @docType package
#' @name RNOmni-help
NULL