# Purpose: Direct INT-based test
# Updated: 19/01/09

#' Direct-INT
#' 
#' Applies the rank-based inverse normal transformation (\code{\link{rankNorm}})
#' to the phenotype \code{y}. Conducts tests of association between the loci in
#' \code{G} and transformed phenotype, adjusting for the model matrix \code{X}.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association. 
#' @param k Offset applied during rank-normalization. See 
#'   \code{\link{rankNorm}}.
#' @param test Either Score or Wald. 
#' @param simple Return the p-values only? 
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @return If \code{simple=T}, returns a vector of p-values, one for each column
#'   of \code{G}. If \code{simple=F}, returns a numeric matrix, including the
#'   Wald or Score statistic, its standard error, the Z-score, and the p-value.
#'   
#' @seealso Basic association test \code{\link{BAT}}, indirect INT \code{\link{IINT}}, omnibus INT \code{\link{OINT}}.
#'   
#' @examples
#' \dontrun{
#' set.seed(100);
#' # Design matrix
#' X = cbind(1,rnorm(1e3));
#' # Genotypes
#' G = replicate(1e3,rbinom(n=1e3,size=2,prob=0.25));
#' storage.mode(G) = "numeric";
#' # Phenotype
#' y = exp(as.numeric(X%*%c(1,1))+rnorm(1e3));
#' # Association test
#' p = DINT(y=y,G=G,X=X,simple=T);
#' }

DINT = function(y,G,X=NULL,k=3/8,test="Score",simple=FALSE,parallel=FALSE){
  # Input check 
  n = length(y);
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(is.null(X)){X=array(1,dim=c(n,1))};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Test
  if(!(test%in%c("Score","Wald"))){stop("Select test from among: Score, Wald.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Please exclude observations missing phenotype or covariate information.")}
  
  # Transform phenotype
  z = rankNorm(y,k=k);
  # Apply basic association test to transformed phenotype
  Out = BAT(y=z,G=G,X=X,test=test,simple=simple);
  return(Out);
}