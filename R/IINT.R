# Purpose: Indirect INT-based method
# Updated: 19/01/11

#' Indirect-INT
#' 
#' Two-stage association testing procedure. In the first stage, phenotype 
#' \code{y} and genotype \code{G} are each regressed on the model matrix
#' \code{X} to obtain residuals. The phenotypic residuals are transformed
#' using \code{\link{rankNorm}}. In the next stage, the INT-transformed
#' residuals are regressed on the genotypic residuals. 
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
#' @param simple Return the p-values only? 
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @return If \code{simple=T}, returns a vector of p-values, one for each column
#'   of \code{G}. If \code{simple=F}, returns a numeric matrix, including the
#'   Wald or Score statistic, its standard error, the Z-score, and the p-value.
#'   
#' @seealso Basic association test \code{\link{BAT}}, direct INT \code{\link{DINT}}, omnibus INT \code{\link{OINT}}.
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
#' p = IINT(y=y,G=G,X=X,simple=T);
#' };

IINT = function(y,G,X=NULL,k=3/8,simple=FALSE,parallel=FALSE){
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
    # SE
    se = sqrt(V);
    # Wald statistic
    U = as.numeric(matIP(g1,e0));
    # Z statistic
    Z = U/se;
    # Chi statistic
    Tw = Z^2;
    # p-value
    p = pchisq(q=Tw,df=1,lower.tail=F);
    # Output
    if(simple){
      Out = c(p);
    } else {
      Out = c(U,se,Z,p);
    }
    return(Out);
  }
  
  # Score statistics
  Out = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  
  ## Format output
  dimnames(Out) = NULL;
  # Check for genotype names
  gnames = colnames(G);
  if(is.null(gnames)){
    gnames = seq(1:ng);
  }
  
  # If returning p-values only
  if(simple){
    names(Out) = gnames;
  } else {
    # Format as matrix
    if(ng==1){
      Out = matrix(Out,nrow=1);
    }
    # Column names
    colnames(Out) = c("Score","SE","Z","p");
    # Row names
    rownames(Out) = gnames;
  };
  # Return
  return(Out);
};