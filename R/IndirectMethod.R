# Purpose: Indirect INT-based method
# Updated: 180516

#' Indirect-INT
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates and structure adjustments to obtain residuals. Genotype is also regressed
#' on covariates and structure adjustments to obtain residuals. In the second stage,
#' INT-transformed phenotypic residuals are regressed on genotypic residuals. 
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq var
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param S Model matrix of structure adjustments.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @param check Logical indicating whether to check the input.
#' @return A numeric matrix of Wald statistics and p-values, one for each locus
#'   (column) in \code{G}, assessing the null hypothesis that genotype is
#'   unrelated to the phenotype.
#'   
#' @examples
#' # IINT against normal phenotype 
#' p = RNOmni::IINT(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

IINT = function(y,G,X,S,k=3/8,parallel=F,check=T){
  if(check){
    # Check inputs
    Input = inCheck(y,G,X,S);
    if(Input$fail){stop("Input check failed.")};
    y = Input$y;
    G = Input$G;
    X = Input$X;
    S = Input$S;
  }
  # Model matrix
  Z = cbind(X,S); 
  # Loci
  n.g = ncol(G);
  # Degrees of freedom
  df2 = length(y)-ncol(Z);
  # Fit null model
  M0 = fitNorm(y=y,Z=Z);
  # Transformed Residuals
  e = rankNormal(M0$Resid);
  # Calculate F statistic
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    e.obs = e[keep];
    # Projected genotype
    Z.obs = Z[keep,,drop=F];
    g.obs = Resid(X=Z.obs,Y=g.obs);
    # Information component
    g2 = sum(g.obs^2);
    # Wald statistic
    r2 = as.numeric(fastIP(A=g.obs,B=e.obs))^2;
    Tw = r2/(g2);
    return(Tw);
  }
  # Wald statistics
  W = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(W,nrow=n.g);
  colnames(Out) = "Wald";
  # Calculate p values
  P = pchisq(q=W,df=1,lower.tail=F);
  Out = cbind(Out,P);
  return(Out);
};