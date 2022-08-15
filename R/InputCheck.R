# Purpose: Inputs checks.
# Updated: 2022-08-15


#' Basic Input Checks
#'
#' Stops evaluation if inputs are improperly formatted.
#' 
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Covariate matrix.
#' @return None.

BasicInputChecks <- function(y, G, X) {
  
  # Ensure y is a numeric vector.
  if (!is.vector(y)) {
    stop("A numeric vector is expected for y.")
  }
  
  # Ensure G is a numeric matrix.
  if (!is.matrix(G)) {
    stop("A numeric matrix is expected for G.")
  }
  
  # Ensure X is a numeric matrix.
  if (!is.matrix(X)) {
    stop("A numeric matrix is expected for X.")
  }
  
  # Ensure y and X are complete.
  is_y_or_x_miss <- any(any(is.na(y)), any(is.na(X)))
  if (is_y_or_x_miss) {
    stop("Please exclude observations missing phenotype or covariate information.")
  }
}