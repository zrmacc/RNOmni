# Purpose: Indirect INT-based method
# Updated: 180912

#' Indirect-INT
#' 
#' Two-stage association testing procedure. In the first stage, phenotype 
#' \code{y} and genotype \code{G} are each regressed on the model matrix
#' \code{X} to obtain residuals. Next, the rank-based inverse normal
#' transformation \code{\link{rankNorm}} is applied to the phenotypic residuals.
#' In the second stage, tests of association are conducted between the genotypic
#' residuals and the transformed phenotypic residuals.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association. 
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNorm}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @return A numeric matrix of Wald statistics and p-values, one for each locus
#'   (column) in \code{G}, assessing the null hypothesis of no genetic effect. 
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
#' p = IINT(y=y,G=G,X=X);
#' };

IINT = function(y,G,X=NULL,k=3/8,parallel=F){
  # Input check 
  n = length(y);
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(is.null(X)){X=array(1,dim=c(n,1))};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Please exclude observations missing phenotype or covariate information.")}
  
  # Loci
  ng = ncol(G);
  # Sample size
  n = length(y);
  # Degrees of freedom
  df2 = n-ncol(X);
  # Fit null model
  M0 = fitOLS(y=y,X=X);
  # Transformed Residuals
  e = rankNorm(M0$Resid);
  # Calculate F statistic
  aux = function(g){
    # Adjust for missingness
    key = !is.na(g);
    miss = (sum(!key)>0);
    if(miss){
      g0 = g[key];
      X0 = X[key,,drop=F];
      e0 = e[key];
    } else {
      g0 = g;
      X0 = X;
      e0 = e;
    }
    # Regression genotype on covariates
    g1 = fitOLS(y=g0,X=X0)$Resid;
    # Information component
    V = sum(g1^2);
    # Wald statistic
    U = as.numeric(matIP(g1,e0));
    Tw = (U^2)/V;
    return(Tw);
  }
  # Wald statistics
  W = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(W,nrow=ng);
  colnames(Out) = "Wald";
  # Calculate p values
  P = pchisq(q=W,df=1,lower.tail=F);
  Out = cbind(Out,P);
  return(Out);
};