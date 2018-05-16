# Purpose: Direct INT-based test
# Updated: 180516

#' Direct-INT
#' 
#' Tests of association between genotype and the rank normalized phenotype, 
#' adjusting for covariates and population structure.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf
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
#' @param check Logical indicating whether to check input formatting.
#' @return A numeric matrix of score statistics and p-values, one for each locus
#'   (column) in \code{G}, assessing the null hypothesis that genotype is
#'   unrelated to the phenotype.
#'   
#' @examples
#' # Direct INT on the normal phenotype 
#' p = RNOmni::DINT(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

DINT = function(y,G,X,S,k=3/8,parallel=F,check=T){
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
  # Sample size
  n = length(y);
  # Degrees of freedom
  df2 = n-ncol(Z);
  # Transform phenotype
  y = rankNormal(y,k=k);
  # Fit null model
  M0 = fitNorm(y=y,Z=Z);
  # Extract model components
  e = M0$Resid;
  tau = M0$Tau;
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    Z.obs = Z[keep,,drop=F];
    e.obs = e[keep];
    # Information components
    I11 = fastIP(g.obs,g.obs);
    I12 = fastIP(g.obs,Z.obs);
    I22 = fastIP(Z.obs,Z.obs);
    # Variance
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    # Score
    a = as.numeric(fastIP(g.obs,e.obs));
    # Test statistic
    Ts = a^2/(V*tau);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(U,nrow=n.g);
  colnames(Out) = "Score";
  # Calculate p values
  P = pf(q=U,df1=1,df2=df2,lower.tail=F);
  Out = cbind(Out,P);
  return(Out);
}