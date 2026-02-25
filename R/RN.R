# Purpose: Rank normal transform
# Updated: 19/10/11

#' Rank-Normalize
#'
#' Applies the rank-based inverse normal transform (INT) to a numeric vector.
#' Observations are first mapped to the (0, 1) scale via the empirical
#' cumulative distribution function (ECDF), then to the real line via the
#' probit (inverse normal CDF).
#'
#' @param u Numeric vector. Must not contain \code{NA}.
#' @param k Offset in the probability scale; must be in (0, 0.5). Default
#'   \code{0.375} corresponds to the Blom transform.
#' @param ties.method Method for breaking ties, passed to \code{\link[base]{rank}}.
#' @return Numeric vector of rank-normalized values (same length as \code{u}).
#' @export
#' @seealso 
#' \itemize{
#'   \item Direct INT test \code{\link{DINT}}.
#'   \item Indirect INT test \code{\link{IINT}}.
#'   \item Omnibus INT test \code{\link{OINT}}.
#' }
#' @examples
#' # Draw from chi-1 distribution
#' y <- stats::rchisq(n = 1e3, df = 1)
#' # Rank normalize
#' z <- RankNorm(y)
#' # Plot density of transformed measurement
#' plot(stats::density(z))
RankNorm <- function(
    u,
    k = 0.375,
    ties.method = "average"
) {
  # Input checks. 
  if (!is.vector(u)) {
    stop("A numeric vector is expected for u.")
  }
  if ((k < 0) || (k > 0.5)) {
    stop("Select the offset within the interval (0,0.5).")
  }
  if (sum(is.na(u)) > 0) {
    stop("Please exclude observations with missing measurements.")
  }

  # Observations.
  n <- length(u)
  
  # Ranks.
  r <- rank(u, ties.method = ties.method)
  
  # Apply transformation.
  out <- stats::qnorm((r - k) / (n - 2 * k + 1))
  return(out)
}
