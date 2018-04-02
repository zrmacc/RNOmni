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
#' @param G Numeric snp by obs genotype matrix.
#' @param X Numeric obs by feature covariate matrix.
#' @param S Numeric obs by feature structure matrix.
#' @param calcP Logical indicating that p values should be calculated.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @return A numeric matrix of score statistics, one for each locus (row) in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{p=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # BAT against normal phenotype
#' p = RNOmni::BAT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

BAT = function(y,G,X,S,calcP=T,parallel=F){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Obs
  n = length(y);
  # Snps
  n.g = nrow(G);
  # Design matrix
  D = cbind(X,S);
  # Degrees of freedom
  df2 = n-ncol(D);
  # Calculate projection matrix
  Q = errProj(D);
  # Calculate scale parameter
  s2 = scaleParam(y=y,Q=Q,df=df2);
  # Caculate score statistic
  aux = function(g){scoreStat(y=y,Q=Q,as.numeric(g),s2=s2)};
  # Calculate score statistics
  U = aaply(.data=G,.margins=1,.fun=aux,.parallel=parallel);
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

#' Direct-INT
#' 
#' Tests for association between genotype and the rank normalized phenotype, 
#' adjusting for covariates and population structure.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Numeric Snp by obs genotype matrix.
#' @param X Numeric Obs by feature covariate matrix.
#' @param S Numeric Obs by feature structure matrix.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See 
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.  
#' @return A numeric matrix of score statistics, one for each locus (row) in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{p=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # Direct INT on the normal phenotype 
#' p = RNOmni::DINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

DINT = function(y,G,X,S,calcP=T,k=3/8,parallel=F){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Obs
  n = length(y);
  # Snps
  n.g = nrow(G);
  # Design matrix
  D = cbind(X,S);
  # Degrees of freedom
  df2 = n-ncol(D);
  q = ncol(D);
  # Transform phenotype
  y = rankNormal(y,k=k);
  # Calculate projection matrix
  Q = errProj(D);
  # Calculate scale parameter
  s2 = scaleParam(y=y,Q=Q,df=df2);
  # Caculate score statistic
  aux = function(g){scoreStat(y=y,Q=Q,g=g,s2=s2)};
  # Calculate score statistics
  U = aaply(.data=G,.margins=1,.fun=aux,.parallel=parallel);
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

#' Fully Indirect-INT
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates and population structure to obtain residuals. In the second stage,
#' INT-transformed residuals are regressed on genotype.
#' 
#' Note that, in simulations, the fully indirect approach did not consistently 
#' control the type I error. For a similar approach that did provide valid 
#' inference, see \code{\link{IINT}}.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf var
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Numeric Snp by obs genotype matrix.
#' @param X Numeric Obs by feature covariate matrix.
#' @param S Numeric Obs by feature structure matrix.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @return A numeric matrix of Wald statistics, one for each locus (row) in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{p=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # FIINT against normal phenotype 
#' p = RNOmni::FIINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

FIINT = function(y,G,X,S,calcP=T,k=3/8,parallel=F){
  warning("This function was included for simulation purposes only.\n");
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Obs
  n = length(y);
  # Snps
  n.g = nrow(G);
  # Design matrix
  D = cbind(X,S);
  # Degrees of freedom
  df2 = n-ncol(D);
  q = ncol(D);
  # Calculate projection matrix
  Q = errProj(D);
  # Transformed Residuals
  e = rankNormal(eps(y=y,Q=Q),k=k);
  # Estimate of scale parameter
  t2 = var(e);
  # Normalize genotypes
  # Calculate F statistic
  aux = function(g){
    # Wald statistic
    g2 = dotP(a=g,b=g);
    r2 = dotP(a=e,b=g)^2;
    Tw = r2/(g2*t2);
    return(Tw);
  }
  # Wald statistics
  W = aaply(.data=G,.margins=1,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(W,nrow=n.g);
  colnames(Out) = "Wald";
  # Calculate p values
  if(calcP){
    P = pf(q=W,df1=1,df2=df2,lower.tail=F);
    Out = cbind(Out,P);
  }
  return(Out);
};

#' Indirect-INT
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates to obtain residuals. In the second stage, INT-transformed
#' residuals are regressed on genotype and population structure.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pf
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @return A numeric vector of p-values assessing the null hypothesis of no 
#'   genotypic effect. P-values are estimated using the Wald statistic, and 
#'   correspond to the rows of G. 
#'   
#' @examples 
#' # IINT against normal phenotype
#' p = RNOmni::IINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

IINT = function(y,G,X,S,calcP=T,k=3/8,parallel=F){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Obs
  n = length(y);
  # Snps
  n.g = nrow(G);
  # Degrees of freedom
  p = ncol(X);
  q = ncol(S);
  df2 = n-p-q;
  # X projection
  Qx = errProj(X);
  # Stage 1 residuals
  ex = rankNormal(eps(y=y,Q=Qx),k=k);
  # S projection 
  Qs = errProj(S);
  # Scale parameter
  s2 = scaleParam(y=ex,Q=Qs,df=df2);
  # Caculate score statistic
  aux = function(g){scoreStat(y=ex,Q=Qs,as.numeric(g),s2=s2)};
  # Calculate score statistics
  U = aaply(.data=G,.margins=1,.fun=aux,.parallel=parallel);
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
  # Estimated correlation
  r = cor(z1,z2);
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
    z2 = -qnorm(IINT(y=yb,G=G,X=Xb,S=Sb,k=k,parallel=F)[,2]);
    return(cbind(z1,z2));
  }
  # Calculate correlation
  aux = function(X){cor(X[1,],X[2,])};
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
  if(p.omni==0){p.omni=1e-16};
  # Output
  Out = c(omni,p.omni);
  return(Out)
};

#' Rank-Normal Omnibus Test
#' 
#' Adaptive association test that synthesizes the \code{\link{DINT}} and 
#' \code{\link{IINT}} approaches. First, the direct and indirect INT based 
#' association tests are applied. An omnibus statistic is calculated based on 
#' whichever approach provides more evidence against the null hypothesis of no 
#' genotypic effect. Details of the method are discussed in the vignette.
#' 
#' Assigning a p-value to the omnibus statistic requires estimation of the 
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
#' @importFrom plyr aaply
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Numeric snp by obs genotype matrix.
#' @param X Numeric obs by feature covariate matrix.
#' @param S Numeric obs by feature structure matrix.
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
#' @return A numeric matrix with three columsn and one row per locus, i.e. row, 
#'   in the genotype matrix, and three columns. The columns are p-values 
#'   obtained by D-INT, IINT, and the omnibus test.
#' 
#' @examples
#' y = RNOmni::Y[,1];
#' Gsub = RNOmni::G[1:10,];
#' X = RNOmni::X;
#' S = RNOmni::S;
#' # Omnibus test against normal phenotype using the average correlation method 
#' p = RNOmni::RNOmni(y=y,G=Gsub,X=X,S=S,method="AvgCorr");
#' # Omnibus test against normal phenotype using the bootstrap correlation method
#' p = RNOmni::RNOmni(y=y,G=Gsub,X=X,S=S,method="Bootstrap",B=10);

RNOmni = function(y,G,X,S,method="AvgCorr",k=3/8,B=100,set.rho,keep.rho=F,keep.stats=F,parallel=F){
  ## Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  y = Input$y;
  G = Input$G;
  X = Input$X;
  S = Input$S;
  ## Additional input checks
  # Dimension check
  ng = nrow(G);
  if(ng<10 & method=="AvgCorr"){stop("Average correlation is not applicable to this few loci.")};
  # Method selection
  flag.m = (method %in% c("AvgCorr","Bootstrap","Manual"));
  if(!flag.m){stop("Select 'AvgCorr', 'Bootstrap', or 'Manual' as the method for correlation estimation.")};
  
  ## Association testing
  # Calculate D-INT p-values
  P1 = suppressWarnings(DINT(y=y,G=G,X=X,S=S,k=k,parallel=parallel));
  # Calculate PI-INT p-values
  P2 = suppressWarnings(IINT(y=y,G=G,X=X,S=S,k=k,parallel=parallel));
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
