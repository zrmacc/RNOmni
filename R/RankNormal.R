# Purpose: Rank normal transform
# Updated: 180516

#' @useDynLib RNOmni
#' @importFrom Rcpp sourceCpp
NULL

#' Rank-Normalize
#' 
#' Applies the rank based inverse normal transform (INT) to a numeric vector. 
#' INT is best-suited to continuous outcomes. See the vignette for the
#' definition of INT.
#' 
#' @importFrom stats qnorm
#' @export
#' 
#' @param u Numeric vector.
#' @param k Offset. Defaults to (3/8), correspond to the Blom transform.
#' @return Numeric vector of rank normalized measurements.
#'   
#' @examples 
#' # Draw from chi-1 distribution
#' y = rchisq(n=1000,df=1);
#' # Rank normalize
#' z = RNOmni::rankNormal(y);
#' # Plot density of transformed measurement
#' plot(density(z));

rankNormal = function(u,k=3/8){
  # Observations
  n = length(u);
  # Ranks
  r = rank(u);
  # Offset
  k = 3/8;
  # Apply transformation
  Out = (r-k)/(n-2*k+1);
  Out = qnorm(Out);
  return(Out);
};

#' Average Correlation Estimate.
#' 
#' Estimates the correlation between correlated p-values on the Z-score scale. The
#' p-values are supposed to have arisen from different tests of the same hypothesis
#' using the same data. Since an estimate of correlation under the null is of interest,
#' pairs where at least one of the Z scores exceeds the threshold \eqn{\tau} are excluded.
#' 
#' @importFrom stats qnorm
#' @export
#' 
#' @param p1 Fist p-value.
#' @param p2 Second p-value.
#' @param tau Threshold Z score above which the p-value likely corresponds to a true positive.
#' @param a Force correlation estimate to fall in the interval (a,1-a);
#' @return A numeric correlation. 

AvgCorr = function(p1,p2,tau=3,a=1e-3){
  # Convert to z-scores
  z1 = -qnorm(p1);
  z2 = -qnorm(p2);
  # Restrict to probable null loci
  keep = (z1<=3)&(z2<=3);
  z1 = z1[keep];
  z2 = z2[keep];
  # Estimated correlation
  r = vecCor(z1,z2);
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
#' @param S Model matrix of structure adjustments.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param B Bootstrap samples for correlation estimation. 
#' @param parallel Run bootstraps in parallel? Must register parallel backend first.
#' @return Numeric matrix of correlation estimates, one per locus (column) in G. 

BootCorr = function(y,G,X,S,k=3/8,B=100,parallel){
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
    Xb = X[Draw,];
    Sb = S[Draw,];
    # Calculate z-statistics
    z1 = -qnorm(DINT(y=yb,G=G,X=Xb,S=Sb,k=k,parallel=F)[,2]);
    z2 = -qnorm(IINT(y=yb,G=G,X=Xb,S=Sb,k=k,parallel=F)[,2]);
    return(cbind(z1,z2));
  }
  # Calculate correlation
  aux = function(X){
    z1 = X[1,];
    z2 = X[2,];
    keep = (z1<=3)&(z2<=3);
    z1 = z1[keep];
    z2 = z2[keep];
    r = vecCor(z1,z2)
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
#' @export
#' 
#' @param u Test statistic.
#' @param r Correlation.
#' @return Numeric p-value. 

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
#' \code{\link{IINT}} approaches. The first approach directly transforms the
#' phenotype, whereas the second approach forms residuals prior to applying the
#' rank normal transformation (\code{\link{rankNormal}}). In the omnibus test,
#' the direct and indirect tests are separately applied. An omnibus statistic is
#' calculated based on whichever approach provides more evidence against the
#' null hypothesis of no genotypic effect. Details of the method are discussed
#' in the vignette.
#' 
#' Assigning a p-value to the omnibus statistic requires an estimate of the 
#' correlation between the test statistics estimated by \code{DINT} and 
#' \code{IINT}. When many loci are under consideration, a computationally 
#' efficient approach is to take the correlation of the observed test statistics
#' across loci (\code{method="AvgCorr"}). Alternatively, when there are fewer 
#' loci, or locus specific estimates are desired, the correlation may be 
#' estimated using bootstrap (\code{method="Bootstrap"}). When using the 
#' bootstrap approach, consider registering a parallel backend and setting 
#' \code{parallel=T}. To manually provide an estimate of the correlation between
#' the test statistics, set (\code{method="Manual"}) and specify 
#' (\code{set.rho}).
#' 
#' @importFrom stats qnorm
#' @importFrom plyr aaply
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param method Method used to estimate correlation for the omnibus test, 
#'   either "AvgCorr", "Bootstrap", or "Manual".
#' @param k Offset applied during rank-normalization. See 
#'   \code{\link{rankNormal}}.
#' @param B If using \code{method=="Bootstrap"}, number of bootstrap samples for
#'   correlation estimation.
#' @param set.rho If using \code{method=="Manual"}, the fixed value of rho, 
#'   either a single value or a vector of length==nrow(G);
#' @param keep.rho Logical indicating whether to return the correlation 
#'   parameter estimated during omnibus calculation. Defaults to FALSE.
#' @param keep.stats Logical indicating whether to return the interim test 
#'   statistics calculated by DINT and IINT. Defaults to FALSE.
#' @param parallel Logical indicating whether to run in parallel. Must register 
#'   parallel backend first.
#' @return A numeric matrix of p values, three for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{keep.stats=T}, the interim test statistics are retained. If
#'   \code{keep.rho=T}, the estimated correlation between the p values provided
#'   by DINT and IINT is retained.
#' 
#' @examples
#' y = RNOmni::Y[,1];
#' Gsub = RNOmni::G[,1:10];
#' X = RNOmni::X;
#' S = RNOmni::S;
#' # Omnibus test against normal phenotype using the average correlation method 
#' p = RNOmni::RNOmni(y=y,G=Gsub,X=X,S=S,method="AvgCorr");
#' # Omnibus test against normal phenotype using the bootstrap correlation method
#' p = RNOmni::RNOmni(y=y,G=Gsub,X=X,S=S,method="Bootstrap",B=10);

RNOmni = function(y,G,X,S,method="AvgCorr",k=3/8,B=100,set.rho,
                  keep.rho=F,keep.stats=F,parallel=F){
  ## Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Mandatory checks
  ng = ncol(G);
  if(ng<10 & method=="AvgCorr"){stop("Average correlation is not applicable to this few loci.")};
  # Method selection
  flag.m = (method %in% c("AvgCorr","Bootstrap","Manual"));
  if(!flag.m){stop("Select 'AvgCorr', 'Bootstrap', or 'Manual' as the method for correlation estimation.")};
  ## Association testing
  # Calculate D-INT p-values
  P1 = DINT(y=y,G=G,X=X,S=S,k=k,parallel=parallel,check=F);
  # Calculate PI-INT p-values
  P2 = IINT(y=y,G=G,X=X,S=S,k=k,parallel=parallel,check=F);
  # Specify correlation between z(DINT) and z(PIINT)
  if(method=="AvgCorr"){
    # Obtain correlation by averaging across loci
    R = AvgCorr(p1=P1[,2],p2=P2[,2]);
  } else if(method=="Bootstrap") {
    # Obtain bootstrap correlation estimate
    R = BootCorr(y=y,G=G,X=X,S=S,B=B,parallel=parallel);
  } else {
    # Manually specify correlation
    if(missing(set.rho)){stop("Provide a value for rho if using method=='Manual'.")}
    R = set.rho;
  }
  # Omnibus statistic
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
  # Keep stats if requested
  if(keep.stats){
    Out = cbind(P1,P2,Q[,1],POmni);
    colnames(Out) = c("DINT.Stat","DINT.p","IINT.Stat","IINT.p","Omni.Stat","Omni.p");
  } else {
    Out = cbind(P1[,2],P2[,2],POmni);
    colnames(Out) = c("DINT","IINT","Omni");
  }
  # Keep correlation if requested 
  if(keep.rho){
    Out = cbind(Out,"Corr"=R);
  }
  return(Out);
}