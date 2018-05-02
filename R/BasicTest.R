# Purpose: Basic score test
# Updated: 180501

#' Basic Association Test
#' 
#' Tests the association between genotype and the untransformed phenotype, 
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
#' @param calcP Logical indicating that p values should be calculated.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @param check Logical indicating whether to check the input. 
#' @return A numeric matrix of score statistics, one for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{calcP=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # BAT against normal phenotype
#' p = RNOmni::BAT(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

BAT = function(y,G,X,S,calcP=T,parallel=F,check=T){
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
  # Fit null model
  M0 = fitNorm(y=y,Z=Z);
  # Extract model components
  e = M0$Resid;
  tau = M0$Tau;
  I22 = tau*M0$Ibb;
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    Z.obs = Z[keep,,drop=F];
    Z.miss = Z[!keep,,drop=F];
    e.obs = e[keep];
    # Information components
    I11 = fastIP(A=g.obs,B=g.obs);
    I12 = fastIP(A=g.obs,B=Z.obs);
    I22.obs = I22-fastIP(Z.miss,Z.miss);
    # Variance
    V = as.numeric(SchurC(I11=I11,I22=I22.obs,I12=I12));
    # Score
    a = as.numeric(fastIP(A=g.obs,B=e.obs));
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
  if(calcP){
    P = pf(q=U,df1=1,df2=df2,lower.tail=F);
    Out = cbind(Out,P);
  }
  return(Out);
}