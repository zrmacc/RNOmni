# Purpose: Omnibus test
# Updated: 180912

#' Average Correlation Estimate.
#' 
#' Estimates the correlation between correlated p-values on the Z-score scale. The
#' p-values are taken as having arisen from different tests of the same hypothesis
#' using the same data. Since an estimate of correlation under the null is of interest,
#' pairs where at least one of the Z scores exceeds the threshold \eqn{\tau} are excluded.
#' 
#' @importFrom stats qnorm
#' 
#' @param p1 Fist p-value.
#' @param p2 Second p-value.
#' @param tau Threshold Z score above which the p-value likely corresponds to a true positive.
#' @param a Force correlation estimate to fall in the interval (a,1-a);
#' @return Numeric correlation. 

AvgCorr = function(p1,p2,tau=3,a=1e-3){
  # Convert to z-scores
  z1 = -qnorm(p1);
  z2 = -qnorm(p2);
  # Restrict to probable null loci
  keep = (z1<=tau)&(z2<=tau);
  z1 = z1[keep];
  z2 = z2[keep];
  # Estimated correlation
  r = cov(A=z1,B=z2,cor=T);
  r = min(r,1-a);
  r = max(r,a);
  # Output
  return(r);
}

#' Bootstrap Correlation Estimate.
#' 
#' Estimates the correlation between correlated p-values on the Z-score scale.
#' Avoids the assumption that the correlation between p-values is constant across
#' loci. Instead, bootstrap is used to calculate locus-specific estimates of the
#' correlation between p-values. 
#' 
#' @importFrom abind abind
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom plyr aaply
#' @importFrom stats cor qnorm
#'  
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param k Offset applied during rank-normalization.
#' @param B Bootstrap samples for correlation estimation. 
#' @param parallel Run bootstraps in parallel? Must register parallel backend first.
#' @return Numeric matrix of correlation estimates, one per locus (column) in G. 

BootCorr = function(y,G,X,k=3/8,B=100,parallel=F){
  # Parallelize
  if(!parallel){foreach::registerDoSEQ()};
  # Obs
  n = length(y);
  # Bind results into an array
  bind3 = function(...){abind::abind(...,along=3)};
  Rho = foreach(i=1:B,.combine="bind3",.multicombine=T,.inorder=F) %dopar% {
    # Draw sample
    Draw = sample(x=n,replace=T);
    # Bootstrap frames
    yb = y[Draw];
    Xb = X[Draw,,drop=F];
    # Calculate z-statistics
    z1 = -qnorm(DINT(y=yb,G=G,X=X,k=k,parallel=F)[,2]);
    z2 = -qnorm(IINT(y=yb,G=G,X=X,k=k,parallel=F)[,2]);
    return(cbind(z1,z2));
  }
  # Calculate correlation
  aux = function(X){
    z1 = X[1,];
    z2 = X[2,];
    keep = (z1<=3)&(z2<=3);
    z1 = z1[keep];
    z2 = z2[keep];
    r = cov(z1,z2,cor=T);
    return(r);
    };
  R = aaply(.data=Rho,.margins=1,.fun=aux,.parallel=parallel);
  # Output
  return(R);
}

#' Omnibus P-value
#' 
#' Calculates the p-value for the maximum of two correlated, standard normal 
#' random variables. 
#' 
#' @importFrom stats dnorm integrate pnorm
#' 
#' @param u Test statistic.
#' @param r Correlation.
#' @return Scalar. 

# Calculate p for omnibus statistic
OmniP = function(u,r){
  # Density of maximum
  f = function(u){
    w = sqrt((1-r)/(1+r));
    return(2*dnorm(u)*pnorm(w*u));
  }
  # P-value
  p = integrate(f=f,lower=u,upper=Inf)$value;
  return(p);
};

#' Rank-Normal Omnibus Test
#' 
#' Association test that synthesizes the \code{\link{DINT}} and 
#' \code{\link{IINT}} methods. The first approach directly transforms the 
#' phenotype, whereas the second approach forms residuals prior to applying the 
#' rank normal transformation (\code{\link{rankNorm}}). In the omnibus test, the
#' direct and indirect tests are separately applied. An omnibus statistic is 
#' calculated based on whichever approach provides more evidence against the 
#' null hypothesis. Details of the method are discussed in the vignette.
#' 
#' Assigning a p-value to the omnibus statistic requires an estimate of the 
#' correlation between the test statistics estimated by \code{DINT} and 
#' \code{IINT} under the null. When many loci are under consideration, a 
#' computationally efficient approach is to take the correlation of the observed
#' test statistics across loci (\code{method="AvgCorr"}). Alternatively, when 
#' there are fewer loci, or locus-specific estimates are desired, the 
#' correlation may be estimated using bootstrap (\code{method="Bootstrap"}). 
#' When using the bootstrap approach, consider registering a parallel backend 
#' and setting \code{parallel=T}. To manually provide an estimate of the 
#' correlation between the test statistics, set (\code{method="Manual"}) and 
#' specify (\code{set.rho}).
#' 
#' @importFrom stats qnorm
#' @importFrom plyr aaply
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by locus genotype matrix.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param method Method used to estimate correlation for the omnibus test, 
#'   either "AvgCorr", "Bootstrap", or "Manual".
#' @param k Offset applied during rank-normalization. See 
#'   \code{\link{rankNorm}}.
#' @param B If using \code{method=="Bootstrap"}, number of bootstrap samples for
#'   correlation estimation.
#' @param set.rho If using \code{method=="Manual"}, the fixed value of rho, 
#'   either a single value or a vector with one element per column in \code{G}.
#' @param keep.rho Logical indicating whether to return the correlation 
#'   parameter estimated during omnibus calculation. Defaults to FALSE.
#' @param keep.stats Logical indicating whether to return the interim test 
#'   statistics calculated by DINT and IINT. Defaults to FALSE.
#' @param parallel Logical indicating whether to run in parallel. Must register 
#'   parallel backend first.
#' @return A numeric matrix of p-values, three for each locus in \code{G}, 
#'   assessing the null hypothesis of no genetic effect. If \code{keep.stats=T},
#'   the interim test statistics are retained. If \code{keep.rho=T}, the
#'   correlation between the p-values provided by DINT and IINT is retained.
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
#' # Average correlation
#' p = RNOmni(y=y,G=G,X=X,method="AvgCorr");
#' # Bootstrap correlation
#' p = RNOmni(y=y,G=G[,1:10],X=X,method="Bootstrap",B=100);
#' # Manual correlation
#' p = RNOmni(y=y,G=G,X=X,method="Manual",set.rho=0.5);
#' }

RNOmni = function(y,G,X=NULL,method="AvgCorr",k=3/8,B=100,set.rho=NULL,keep.rho=F,keep.stats=F,parallel=F){
  ## Check inputs
  # Input check 
  n = length(y);
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  if(is.null(X)){X=array(1,dim=c(n,1))};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Please exclude observations missing phenotype or covariate information.")}
  
  # Correlation calculation
  ng = ncol(G);
  if((ng<25)&&(method=="AvgCorr")){stop("Average correlation is not applicable to this few loci.")};
  # Method selection
  Choices = c("AvgCorr","Bootstrap","Manual");
  if(!(method%in%Choices)){stop("Select 'AvgCorr', 'Bootstrap', or 'Manual' as the method for correlation estimation.")};
  if((method=="Manual")){
    if(is.null(set.rho)||(min(set.rho)<=0)||(max(set.rho)>=1)){
      stop("If using 'Manual' correlction, initialize the element(s) of set.rho in the interval (0,1).")};
  }
  
  ## Association testing
  # Calculate D-INT p-values
  P1 = DINT(y=y,G=G,X=X,k=k,parallel=parallel);
  # Calculate PI-INT p-values
  P2 = IINT(y=y,G=G,X=X,k=k,parallel=parallel);
  # Specify correlation between z(DINT) and z(PIINT)
  if(method=="AvgCorr"){
    # Obtain correlation by averaging across loci
    R = AvgCorr(p1=P1[,2],p2=P2[,2]);
  } else if(method=="Bootstrap") {
    # Obtain bootstrap correlation estimate
    R = BootCorr(y=y,G=G,X=X,B=B,parallel=parallel);
  } else {
    # Manually specify correlation
    R = set.rho;
  }
  
  ## Omnibus statistic
  Q = cbind(-qnorm(P1[,2]),-qnorm(P2[,2]));
  Q = aaply(.data=Q,.margins=1,.fun=max);
  Q = cbind(Q,R);
  # Calculate Omnibus p-values
  aux = function(v){
    OmniP(v[1],v[2]);
  }
  if(ng==1){
    POmni = matrix(aux(Q),nrow=1);
  } else {
    POmni = aaply(.data=Q,.margins=1,.fun=aux,.parallel=parallel);
  }
  
  ## Keep stats if requested
  if(keep.stats){
    Out = cbind(P1,P2,Q[,1],POmni);
    colnames(Out) = c("DINT.Stat","DINT.p","IINT.Stat","IINT.p","Omni.Stat","Omni.p");
  } else {
    Out = cbind(P1[,2],P2[,2],POmni);
    colnames(Out) = c("DINT","IINT","Omni");
  }
  
  ## Keep correlation if requested 
  if(keep.rho){
    Out = cbind(Out,"Corr"=R);
  }
  return(Out);
}