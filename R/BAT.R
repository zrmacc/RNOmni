# Purpose: Basic score test
# Updated: 180912

#' Basic Association Test
#' 
#' Conducts tests of association between the loci in \code{G} and the
#' untransformed phenotype \code{y}, adjusting for the model matrix \code{X}.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association. 
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @return A numeric matrix of score statistics and p-values, one for each locus
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
#' y = as.numeric(X%*%c(1,1))+rnorm(1e3);
#' # Association test
#' p = BAT(y=y,G=G,X=X);
#' }

BAT = function(y,G,X=NULL,parallel=F){
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
  # Degrees of freedom
  df2 = n-ncol(X);
  # Fit null model
  M0 = fitOLS(y=y,X=X);
  # Extract model components
  e = M0$Resid;
  v = M0$V;
  # Base information
  Iaa = M0$Ibb*v;
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    key = !is.na(g);
    miss = (sum(!key)>0);
    if(miss){
      g0 = g[key];
      X0 = X[key,,drop=F];
      X1 = X[!key,,drop=F];
      e0 = e[key];
    } else {
      g0 = g;
      X0 = X;
      e0 = e;
    }
    # Information components
    I11 = matIP(g0,g0);
    I12 = matIP(g0,X0);
    I22 = Iaa;
    if(miss){
      # Loss of information
      I22 = I22-matIP(X1,X1);
    }
    # Efficient info
    E = as.numeric(SchurC(I11,I22,I12));
    # Score
    U = as.numeric(matIP(g0,e0));
    # Test statistic
    Ts = (U^2)/(v*E);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(U,nrow=ng);
  colnames(Out) = "Score";
  # Calculate p values
  P = pf(q=U,df1=1,df2=df2,lower.tail=F);
  Out = cbind(Out,P);
  return(Out);
}