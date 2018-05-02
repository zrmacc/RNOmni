#' @useDynLib RNOmni
#' @importFrom Rcpp sourceCpp
NULL

#' Rank-Normalize
#' 
#' Applies the rank based inverse normal transform (INT) to a numeric vector. 
#' INT is indicated for continuous outcomes. See the vignette for the
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
#' Estimate correlation using the average of qnorm(p1)*qnorm(p2) across loci, where
#' (p1,p2) are p-values obtained via two different association tests. 
#' 
#' @param p1 Fist p-value.
#' @param p2 Second p-value.
#' @param a Force correlation estimate to fall in the interval (a,1-a);
#' @importFrom stats qnorm

AvgCorr = function(p1,p2,a=1e-3){
  # Convert to z-scores
  z1 = -qnorm(p1);
  z2 = -qnorm(p2);
  # Restrict to probable null loci
  keep = (abs(z1)<=3)&(abs(z2)<=3);
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
#' Use bootstrap to estimate correlation among Z statistics
#' whose maximum is taken in the omnibus test.
#' 
#' @importFrom abind abind
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom plyr aaply
#' @importFrom stats cor qnorm
#'  
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param B Bootstrap samples for correlation estimation. 
#' @param parallel Run bootstraps in parallel? Must register parallel backend first.

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
    z2 = -qnorm(IINTd(y=yb,G=G,X=Xb,S=Sb,k=k,parallel=F)[,2]);
    return(cbind(z1,z2));
  }
  # Calculate correlation
  aux = function(X){
    z1 = X[1,];
    z2 = X[2,];
    keep = (abs(z1)<=3)&(abs(z2)<=3);
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
#' Calculate p-value for omnibus statistic.
#' @param Q Numeric vector formatted as (p1,p2,rho), where p1 and p2 are 
#' estimated p-values, and rho is their correlation. 
#' @importFrom stats qnorm
#' @importFrom mvtnorm pmvnorm
#' 
#' @export

# Calculate p for omnibus statistic
OmniP = function(Q){
  # Q : vector formatted as (p1,p2,rho);
  p = Q[1:2];
  r = Q[3];
  # Convert to Z's
  Z = -qnorm(p);
  # Omnibus statistic
  omni = max(Z);
  # Correlation matrix
  C = matrix(c(1,r,r,1),ncol=2,byrow=T);
  # P-value
  p.omni = 1 - mvtnorm::pmvnorm(upper=c(omni,omni),corr=C)[1];
  if(p.omni<=1e-16){p.omni=1e-16};
  # Output
  Out = c(omni,p.omni);
  return(Out)
};

#' Rank-Normal Omnibus Test
#' 
#' Association test that synthesizes the \code{\link{DINT}} and 
#' \code{\link{IINTc}} approaches. First, the direct and indirect association
#' tests are applied. An omnibus statistic is calculated based on whichever
#' approach provides more evidence against the null hypothesis of no genotypic
#' effect. Details of the method are discussed in the vignette.
#' 
#' Assigning a p-value to the omnibus statistic requires estimation of the 
#' correlation between the test statistics estimated by \code{DINT} and 
#' \code{IINTc}. When many loci are under consideration, a computationally 
#' efficient approach is to take the correlation of the observed test statistics
#' across loci (\code{method="AvgCorr"}). Alternatively, when there are fewer 
#' loci, or locus specific estimates are desired, the correlation may be 
#' estimated using bootstrap (\code{method="Bootstrap"}). When using the 
#' bootstrap approach, consider registering a parallel backend and setting 
#' \code{parallel=T}. To manually provide an estimate of the correlation between
#' the test statistics, set (\code{method="Manual"}) and specify 
#' (\code{set.rho}).
#' 
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
#' @param check Logical indicating whether to check the input.
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
                  keep.rho=F,keep.stats=F,parallel=F,check=T){
  if(check){
    ## Check inputs
    Input = inCheck(y,G,X,S);
    if(Input$fail){stop("Input check failed.")};
    y = Input$y;
    G = Input$G;
    X = Input$X;
    S = Input$S;
  }
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
  P2 = IINTd(y=y,G=G,X=X,S=S,k=k,parallel=parallel,check=F);
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
  # Matrix containing P1, P2, and their estimated correlation;
  Q = cbind("p1"=P1[,2],"p2"=P2[,2],R);
  # Calculate Omnibus p-values
  if(ng==1){
    POmni = matrix(OmniP(Q),nrow=1);
  } else {
    POmni = aaply(.data=Q,.margins=1,.fun=OmniP,.parallel=parallel);
  }
  # Keep stats if requested
  if(keep.stats){
    Out = cbind(P1,P2,POmni);
    colnames(Out) = c("DINT.Score","DINT.p","IINT.Score","IINT.p","Omni.Stat","Omni.p");
  } else {
    Out = cbind(P1[,2],P2[,2],POmni[,2]);
    colnames(Out) = c("DINT","IINT","Omni");
  }
  # Keep correlation if requested 
  if(keep.rho){
    Out = cbind(Out,"Corr"=R);
  }
  return(Out);
}
