# Purpose: Basic score test
# Updated: 181029

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
#' @param test Either Score or Wald. 
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

BAT = function(y,G,X=NULL,test="Score",parallel=F){
  # Input check 
  n = length(y);
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(is.null(X)){X=array(1,dim=c(n,1))};
  k = ncol(X);
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Test
  if(!(test%in%c("Score","Wald"))){stop("Select test from among: Score, Wald.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Please exclude observations missing phenotype or covariate information.")};
  
  # Loci
  ng = ncol(G);

  ## Score test
  if(test=="Score"){
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
      # p-value
      p = pf(q=Ts,df1=1,df2=sum(key)-k-1,lower.tail=F);
      # Output
      Out = c(Ts,p);
      return(Out);
    }
    # Calculate score statistics
    Out = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel,.drop=F);
  } else {
  ## Wald Test
    # Function to calculate wald statistics
    aux = function(g){
      # Adjust for missingness
      key = !is.na(g);
      miss = (sum(!key)>0);
      if(miss){
        y0 = y[key];
        g0 = g[key];
        X0 = X[key,,drop=F];
      } else {
        y0 = y;
        g0 = g;
        X0 = X;
      }
      # Fit full model
      M1 = fitOLS(y=y0,X=cbind(g0,X0));
      # Coefficient
      bg = M1$Beta[1];
      # Variance
      Ibbi = as.numeric(matInv(M1$Ibb)[1,1]);
      # Test statistic
      Ts = (bg^2)/(Ibbi);
      # p-value
      p = pf(q=Ts,df1=1,df2=sum(key)-k-1,lower.tail=F);
      # Output
      Out = c(Ts,p);
      return(Out);
    }
    # Calculate wald statistics
    Out = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel,.drop=F);
  }
  # Format
  dimnames(Out) = NULL;
  colnames(Out) = c(test,"p");
  if(!is.null(colnames(G))){rownames(Out)=colnames(G)} else {rownames(Out) = seq(1:ng)};
  # Return
  return(Out);
}