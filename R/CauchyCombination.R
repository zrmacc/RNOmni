# Purpose: Cauchy combination.
# Updated: 2022-08-15


#' Convert P-value to Cauchy Random
#' 
#' @param p Numeric p-value.
#' @return Numeric Cauchy random variable. 
PtoCauchy <- function(p) {
  if (p < 1e-10) {
    out <- 1 / (p * pi)
  } else {
    out <- tanpi(0.5 - p)
  }
  return(out)
}


#' Convert Cauchy Random Variable to P
#' 
#' @param z Numeric Cauchy random variable. 
#' @return Numeric p-value.

CauchyToP <- function(z) {
  if (z > 1e10) {
    out <- (1 / z) / pi
  } else {
    out <- stats::pcauchy(q = z, lower.tail = FALSE)
  }
  return(out)
}


# -----------------------------------------------------------------------------

#' Omnibus P-value.
#' 
#' Obtains an omnibus p-value from a vector of potentially dependent p-values 
#' using the method of Cauchy combination. The p-values are converted to Cauchy
#' random deviates then averaged. The distribution of the average of these 
#' deviates is well-approximated by a Cauchy distribution in the tails. See
#' <https://doi.org/10.1080/01621459.2018.1554485>.
#'
#' @param p Numeric vector of p-values.
#' @param w Numeric weight vector.
#' @return OINT p-value.
#' @export

OmniP <- function(p, w = NULL) {
  # Check input
  m <- min(p)
  M <- max(p)
  
  if (is.null(w)) {
    w <- rep(1, length(p))
  } else {
    are_equal <- all.equal(length(p), length(w))
    if (!are_equal) {
      stop("w and p must have the same length.")
    }
  }
  
  # Cases
  if ((m < 0) | (M > 1)) {
    stop("Cannot have p-values < 0 or > 1.")
  }
  if ((m == 0) & (M == 1)) {
    stop("Cannot have p-values of 0 and 1 simultaneously.")
  }
  if ((m == 0) & (M < 1)) {
    return(0)
  }
  if ((m > 0) & (M == 1)) {
    return(1)
  }
  
  # Convert to Cauchy.
  q <- sapply(p, PtoCauchy)
  
  # Mean of Cauchy random variables, which is again Cauchy.
  z <- sum(w * q)
  
  # Invert to obtain p-value.
  p <- CauchyToP(z / sum(w))
  return(p)
}
