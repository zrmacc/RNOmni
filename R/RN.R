# Purpose: Rank normal transform
# Updated: 19/10/11

#' Rank-Normalize
#' 
#' Applies the rank based inverse normal transform (INT) to a numeric vector.
#' The INT can be broken down into a two-step procedure. In the first, the
#' observations are transformed onto the probability scale using the empirical
#' cumulative distribution function (ECDF). In the second, the observations are
#' transformed onto the real line, as Z-scores, using the probit function.
#' 
#' @importFrom stats qnorm
#' @export
#' 
#' @param u Numeric vector.
#' @param k Offset. Defaults to (3/8), correspond to the Blom transform.
#' @return Numeric vector of rank normalized measurements.
#' 
#' @seealso Direct INT \code{\link{DINT}}, indirect INT \code{\link{IINT}}, omnibus INT \code{\link{OINT}}.
#'   
#' @examples 
#' \dontrun{
#' # Draw from chi-1 distribution
#' y = rchisq(n=1e3,df=1);
#' # Rank normalize
#' z = rankNorm(y);
#' # Plot density of transformed measurement
#' plot(density(z));
#' }

rankNorm = function(u,k=3/8){
  if(!is.vector(u)){stop("A numeric vector is expected for u.")};
  if((k<0)||(k>0.5)){stop("Select the offset within the interval (0,0.5).")};
  if(sum(is.na(u))>0){stop("Please exclude observations with missing measurements.")};
  
  # Observations
  n = length(u);
  # Ranks
  r = rank(u);
  # Apply transformation
  return(qnorm((r-k)/(n-2*k+1)));
};