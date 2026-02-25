# Purpose: Rank Normal Omnibus test.
# Updated: 2022-08-15

# -----------------------------------------------------------------------------

#' Omnibus-INT
#'
#' Association test that synthesizes the \code{\link{DINT}} and
#' \code{\link{IINT}} tests. The first approach is most powerful for traits that
#' could have arisen from a rank-preserving transformation of a latent normal
#' trait. The second approach is most powerful for traits that are linear in
#' covariates, yet have skewed or kurtotic residual distributions. During the
#' omnibus test, the direct and indirect tests are separately applied, then the
#' p-values are combined via the Cauchy combination method.
#'
#' @param y Numeric phenotype vector.
#' @param G Genotype matrix with observations as rows, SNPs as columns.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{RankNorm}}.
#' @param ties.method Method of breaking ties, passed to \code{base::rank}.
#' @param weights Numeric length-2 vector of weights for the D-INT and I-INT
#'   p-values in the Cauchy combination. Default \code{c(1, 1)} gives equal weight.
#' @param simple If \code{TRUE}, return only the omnibus p-values; if
#'   \code{FALSE}, return a matrix with D-INT, I-INT, and OINT p-values.
#' @return If \code{simple = TRUE}, a named numeric vector of omnibus p-values
#'   (one per column of \code{G}). If \code{simple = FALSE}, a numeric matrix
#'   with columns \code{DINT-p}, \code{IINT-p}, \code{OINT-p} and one row per SNP.
#' @export
#' @seealso
#' \itemize{
#'   \item Basic association test \code{\link{BAT}}.
#'   \item Direct INT test \code{\link{DINT}}.
#'   \item Indirect INT test \code{\link{IINT}}.
#' }
#'
#' @examples
#' set.seed(100)
#' # Design matrix
#' X <- cbind(1, rnorm(1e3))
#' # Genotypes
#' G <- replicate(1e3, rbinom(n = 1e3, size = 2, prob = 0.25))
#' storage.mode(G) <- "numeric"
#' # Phenotype
#' y <- exp(as.numeric(X %*% c(1, 1)) + rnorm(1e3))
#' # Omnibus
#' p <- OINT(y = y, G = G, X = X, simple = TRUE)
OINT <- function(
    y,
    G,
    X = NULL,
    k = 0.375,
    ties.method = "average",
    weights = c(1, 1),
    simple = FALSE
) {
  
  # Generate X if omitted.
  if (is.null(X)) {
    X <- array(1, dim = c(length(y), 1))
  }
  
  # Input check.
  BasicInputChecks(y, G, X)

  # Association testing.
  # Calculate D-INT p-values.
  p_dint <- DINT(
    y = y, G = G, X = X, k = k, ties.method = ties.method, simple = TRUE)
  
  # Calculate I-INT p-values.
  p_iint <- IINT(
    y = y, G = G, X = X, k = k, ties.method = ties.method, simple = TRUE)

  # P Matrix
  p_mat <- cbind(p_dint, p_iint)

  # Omnibus p-values
  p_omni <- plyr::aaply(
    .data = p_mat,
    .margins = 1,
    .fun = function(x) {OmniP(p = x, w = weights)}
  )

  # Output

  # Check for genotype names.
  gnames <- colnames(G)
  if (is.null(gnames)) {
    gnames <- seq_len(ncol(G))
  }

  # Format
  if (simple) {
    out <- p_omni
    names(out) <- gnames
  } else {
    out <- cbind(
      "DINT-p" = p_dint,
      "IINT-p" = p_iint, 
      "OINT-p" = p_omni
    )
    rownames(out) <- gnames
  }
  return(out)
}
