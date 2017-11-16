#' Rank-Normalize
#' 
#' Applies the rank based inverse normal transform (INT) to a numeric vector.
#' INT is indicated for continuous phenotypes lacking ties. See the vignette for
#' the mathematical definition of INT.
#' 
#' @importFrom stats qnorm
#' @export
#' 
#' @param u Numeric vector.
#' @param c Offset. Defaults to (3/8), correspond to the Blom transform.
#' @return Numeric vector of rank normalized measurements.
#'   
#' @examples 
#' # Draw from chi-1 distribution
#' y = rchisq(n=1000,df=1);
#' # Rank normalize
#' z = RNOmni::rankNormal(y);
#' # Plot density of transformed measurement
#' plot(density(z));

rankNormal = function(u,c=3/8){
  # Observations
  n = length(u);
  # Ranks
  r = rank(u);
  # Offset
  c = 3/8;
  # Apply transformation
  Out = (r-c)/(n-2*c+1);
  Out = qnorm(Out);
  return(Out);
};

#' Missingness Filter
#' 
#' Function to adjust for missing data. Observations with phenotype or structure
#' adjustments missing are removed. Missing covariates are imputed to the median
#' of the observed values. An observation missing genotype information is
#' excluded from association testing only at those loci where genotype is
#' unobserved.
#' 
#' @importFrom stats median
#' 
#' @param y Numeric phenotype vector
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.


missFilter = function(y,G,X,S){
  # Impute missing X
  aux = function(col){
    col[is.na(col)] = median(col[!is.na(col)]);
    return(col);
  }
  X = apply(X,MARGIN=2,FUN=aux);
  # Subjects with y or S missing
  ind.1 = is.na(y);
  ind.2 = apply(S,MARGIN=1,FUN=function(row){sum(is.na(row))>0});
  keep = !(ind.1|ind.2);
  # Filtering
  y.out = y[keep];
  G.out = G[,keep,drop=F];
  X.out = X[keep,];
  S.out = S[keep,];
  if(sum(keep)<length(y)){warning(sprintf("%i observations remain after missingness filter.\n",sum(keep)));}
  # Output
  Out = list("y"=y.out,"G"=G.out,"X"=X.out,"S"=S.out);
  return(Out);
}

#' Input Check
#' 
#' Function to ensure the dimensions of inputs to association methods agree.
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.

inCheck = function(y,G,X,S){
  # Check phenotype
  flag.y = !(is.numeric(y)&is.vector(y));
  if(flag.y){
    warning("A numeric vector is required for y.");
  }
  n.y = length(y);
  # Check genotype
  flag.g = is.vector(G);
  if(flag.g){
    warning("A matrix is expected for G.");
    G = matrix(G,nrow=1);
  }
  n.g = ncol(G);
  # Check covariates
  flag.x = is.vector(X);
  if(flag.x){
    warning("A matrix or data.frame is expected for X.");
    X = data.frame(X);
  }
  n.x = nrow(X);
  # Check structure matrix
  flag.s = is.vector(S);
  if(flag.s){
    warning("A matrix or data.frame is expected for S.");
    S = data.frame(S);
  }
  n.s = nrow(S);
  # Dimensional consistency
  flag.d = !all.equal(n.y,n.g,n.x,n.s);
  if(flag.d){
    warning("Dimensions of inputs are inconsistent. Ensure length(y)=ncol(G)=nrow(X)=nrow(S).")
  }
  # Output
  Out = list("fail"=flag.y|flag.d,"G"=G,"X"=X,"S"=S);
  return(Out);
}

#' Basic Association Test
#' 
#' Regression of the untransformed phenotype on genotype, covariates, and 
#' adjustments for population structure.
#' 
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param M Apply missingness filter? See \code{\link{missFilter}}.
#' @return A numeric vector of p-values assessing the null hypothesis of no 
#'   genotypic effect. P-values are estimated using the Wald statistic, and 
#'   correspond to the rows of G. 
#'   
#' @examples
#' # BAT against normal phenotype
#' p = RNOmni::BAT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

BAT = function(y,G,X,S,M=T){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Missingness filter
  if(M){
    Miss = missFilter(y,G,X,S);
    y = Miss[["y"]];
    G = Miss[["G"]];
    X = Miss[["X"]];
    S = Miss[["S"]];
  }
  # Design matrix
  D = cbind(X,S);
  D = model.matrix(~.,data=data.frame(D));
  q = ncol(D);
  # Function to implement regression and calculate wald p-value
  getP = function(g){
    # Add genotype
    Dg = cbind(D,g);
    # Missing genotype
    keep = !is.na(g);
    y = y[keep];
    Dg = Dg[keep,];
    # Regression
    M = RcppEigen::fastLmPure(X=Dg,y=y);
    # Wald statistic
    w = coef(M)[q+1]/(M$se[q+1]);
    # P-value
    p = 2*pt(q=abs(w),df=M$df.residual,lower.tail=F)
    return(p);
  }
  # Calculate P-values
  P = apply(X=G,MARGIN=1,FUN=getP);
  return(P);
}

#' Direct-INT
#' 
#' Rank-normalizes the phenotype, then regresses the transformed phenotype on
#' genotype, covariates, and adjustments for population structure.
#'
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix. 
#' @param M Apply missingness filter? See \code{\link{missFilter}}.
#' @param c Offset applied during rank-normalization. See \code{\link{rankNormal}}.
#' @return A numeric vector of p-values assessing the null hypothesis of no 
#'   genotypic effect. P-values are estimated using the Wald statistic, and 
#'   correspond to the rows of G. 
#' 
#' @examples
#' # DINT against normal phenotype 
#' p = RNOmni::DINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

DINT = function(y,G,X,S,M=T,c=3/8){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Missingness filter
  if(M){
    Miss = missFilter(y,G,X,S);
    y = Miss[["y"]];
    G = Miss[["G"]];
    X = Miss[["X"]];
    S = Miss[["S"]];
  }
  # Design matrix
  D = cbind(X,S);
  D = model.matrix(~.,data=data.frame(D));
  q = ncol(D);
  # Transform phenotype
  y = rankNormal(y,c=c);
  # Function to implement regression and calculate wald p-value
  getP = function(g){
    # Add genotype
    Dg = cbind(D,g);
    # Missing genotype
    keep = !is.na(g);
    y = y[keep];
    Dg = Dg[keep,];
    # Regression
    M = RcppEigen::fastLmPure(X=Dg,y=y);
    # Wald statistic
    w = coef(M)[q+1]/(M$se[q+1]);
    # P-value
    p = 2*pt(q=abs(w),df=M$df.residual,lower.tail=F)
    return(p);
  }
  # Calculate P-values
  P = apply(X=G,MARGIN=1,FUN=getP);
  return(P);
}

#' Fully Indirect-INT
#' 
#' Two-stage regression procedure In the first stage, phenotype is regressed on 
#' covariates and adjustments for population structure to obtain residuals. In 
#' the second stage, INT-transformed residuals are regressed on genotype.
#' 
#' Note that, in simulations, FIINT did not consistently provide valid
#' inference. For a similar approach that did control the type I error, see
#' \code{\link{PIINT}}.
#' 
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt resid
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param M Apply missingness filter? See \code{\link{missFilter}}.
#' @param c Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @return A numeric vector of p-values assessing the null hypothesis of no 
#'   genotypic effect. P-values are estimated using the Wald statistic, and 
#'   correspond to the rows of G.
#'   
#' @examples
#' # FIINT against normal phenotype 
#' p = RNOmni::FIINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

FIINT = function(y,G,X,S,M=T,c=3/8){
  warning("This function was included for simulation purposes, and was not found to provide valid inference.\n");
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Missingness filter
  if(M){
    Miss = missFilter(y,G,X,S);
    y = Miss[["y"]];
    G = Miss[["G"]];
    X = Miss[["X"]];
    S = Miss[["S"]];
  }
  # Design matrix
  D = cbind(X,S);
  D = model.matrix(~.,data=data.frame(D));
  # Residuals
  e = resid(RcppEigen::fastLmPure(X=D,y=y));
  # Transform
  z = rankNormal(e,c=c);
  # Function to implement regression and calculate wald p-value
  getP = function(g){
    # Add genotype
    Dg = cbind(1,g);
    # Missing genotype
    keep = !is.na(g);
    z = z[keep];
    Dg = Dg[keep,];
    # Regression
    M = RcppEigen::fastLmPure(X=Dg,y=z);
    # Wald statistic
    w = coef(M)[2]/(M$se[2]);
    # P-value
    p = 2*pt(q=abs(w),df=M$df.residual,lower.tail=F)
    return(p);
  }
  # Calculate P-values
  P = apply(X=G,MARGIN=1,FUN=getP);
  return(P);
};

#' Partially Indirect-INT
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates to obtain residuals. In the second stage, INT-transformed
#' residuals are regressed on genotype and adjustments for population structure.
#' 
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt resid
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param M Apply missingness filter? See \code{\link{missFilter}}.
#' @param c Offset applied during rank-normalization. See \code{\link{rankNormal}}.
#' @return A numeric vector of p-values assessing the null hypothesis of no 
#'   genotypic effect. P-values are estimated using the Wald statistic, and 
#'   correspond to the rows of G. 
#'   
#' @examples 
#' # PIINT against normal phenotype
#' p = RNOmni::PIINT(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S);

PIINT = function(y,G,X,S,M=T,c=3/8){
  # Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  G = Input$G;
  X = Input$X;
  S = Input$S;
  # Missingness filter
  if(M){
    Miss = missFilter(y,G,X,S);
    y = Miss[["y"]];
    G = Miss[["G"]];
    X = Miss[["X"]];
    S = Miss[["S"]];
  }
  # Design matrix
  D = model.matrix(~.,data=data.frame(X));
  # Residuals
  e = resid(RcppEigen::fastLmPure(X=D,y=y));
  # Transform
  z = rankNormal(e,c=c);
  # Function to implement regression and calculate wald p-value
  getP = function(g){
    # Add genotype
    Dg = cbind(1,g,S);
    # Missing genotype
    keep = !is.na(g);
    z = z[keep];
    Dg = Dg[keep,];
    # Regression
    M = RcppEigen::fastLmPure(X=Dg,y=z);
    # Wald statistic
    w = coef(M)[2]/(M$se[2]);
    # P-value
    p = 2*pt(q=abs(w),df=M$df.residual,lower.tail=F)
    return(p);
  }
  # Calculate P-values
  P = apply(X=G,MARGIN=1,FUN=getP);
  return(P);
}

#' Average Correlation Estimate.
#' 
#' Estimate correlation using the average of qnorm(p1)*qnorm(p2) across loci, where
#' (p1,p2) are p-values obtained via two different association tests. 
#' 
#' @param P1 Fist p-value.
#' @param P2 Second p-value.
#' @param eps Force correlation estimate to fall in the interval (eps,1-eps);
#' @importFrom stats qnorm

AvgCorr = function(P1,P2,eps=1e-3){
  ng = length(P1);
  # Convert to z-scores
  z1 = -qnorm(P1);
  z2 = -qnorm(P2);
  # Estimated correlation
  r = mean(z1*z2);
  r = min(r,1-eps);
  r = max(r,eps);
  # Output
  R = rep(r,times=ng);
  return(R);
}

#' Bootstrap Correlation Estimate.
#' 
#' Use bootstrap to estimate correlation among Z statistics
#' whose maximum is taken in the omnibus test. 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param B Bootstrap samples for correlation estimation. 
#' @param parallel Run bootstraps in parallel? Must register parallel backend first.
#' @importFrom abind abind
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom stats cor qnorm


BootCorr = function(y,G,X,S,B=100,parallel){
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
    z1 = -qnorm(DINT(y=yb,G=G,X=Xb,S=Sb,M=F));
    z2 = -qnorm(PIINT(y=yb,G=G,X=Xb,S=Sb,M=F));
    return(cbind(z1,z2));
  }
  # Calculate correlation
  aux = function(X){cor(X[1,],X[2,])};
  R = apply(X=Rho,MARGIN=1,FUN=aux);
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

# Calculate p for omnibus statistic
Omni = function(Q){
  # Q : vector formatted as (p1,p2,rho);
  p = Q[1:2];
  r = Q[3];
  # Convert to Z's
  Z = -qnorm(p);
  # Omnibus statistic
  o = max(Z);
  # Correlation matrix
  C = matrix(c(1,r,r,1),ncol=2,byrow=T);
  # P-value
  Out = 1 - mvtnorm::pmvnorm(upper=c(o,o),corr=C)[1];
  return(Out)
};

#' Rank-Normal Omnibus Test
#' 
#' Omnibus association test that synthesizes the \code{\link{DINT}} and
#' \code{\link{PIINT}} approaches. In the omnibus test, both DINT and PIINT are
#' applied. An omnibus statistic is calculated based on whichever approach
#' provides more evidence against the null hypothesis of no genotypic effect.
#' Details of the method are discussed below and in the vignette.
#' 
#' Assignment of a p-value to the omnibus statistic requires an estimate of the
#' correlation between the test statistics estimated by DINT and PIINT. When the
#' sample size and number of loci are both large, and efficient estimate of the
#' correlation is obtained by averaging across loci (\code{method="AvgCorr"}).
#' When either the sample size or the number of loci is small, bootstrap
#' (\code{method="Bootstrap"}) allows for locus specific correlation estimates.
#' If using the bootstrap approach, consider registering a parallel backend and
#' setting \code{parallel=T}. To manually provide an estimate of the correlation
#' between the test statistics, set (\code{method="Manual"}) and specify (\code{set.rho}).
#' 
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param method Method used to estimate correlation for the omnibus test, 
#'   either "AvgCorr", "Bootstrap", or "Manual".
#' @param M Apply missingness filter? See \code{\link{missFilter}}.
#' @param c Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param B If using \code{method=="Bootstrap"}, number of bootstrap samples for
#'   correlation estimation.
#' @param set.rho If using \code{method=="Manual"}, the fixed value of rho,
#'   either a single value or a vector of length==nrow(G);
#' @param keep.rho Logical indicating whether to return the correlation parameter 
#'   estimated during omnibus calculation. Defaults to FALSE.
#' @param parallel Run bootstraps in parallel? Must register parallel backend
#' first.
#' @return A numeric matrix with three columsn and one row per locus, i.e. row,
#'   in the genotype matrix, and three columns. The columns are p-values
#'   obtained by DINT, PIINT, and the omnibus test.
#'   
#' @examples
#' # Omnibus test against normal phenotype using the average correlation method 
#' p = RNOmni::RNOmni(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S,method="AvgCorr");
#' # Omnibus test against normal phenotype using the bootstrap correlation method
#' p = RNOmni::RNOmni(y=RNOmni::Y[,1],G=RNOmni::G[1:10,],X=RNOmni::X,S=RNOmni::S,method="Bootstrap");

RNOmni = function(y,G,X,S,method="AvgCorr",M=T,c=3/8,B=100,set.rho,keep.rho=F,parallel=F){
  ## Check inputs
  Input = inCheck(y,G,X,S);
  if(Input$fail){stop("Input check failed.")};
  G = Input$G;
  X = Input$X;
  S = Input$S;
  
  ## Missingness Filter 
   if(M){
    Miss = missFilter(y,G,X,S);
    y = Miss[["y"]];
    G = Miss[["G"]];
    X = Miss[["X"]];
    S = Miss[["S"]];
  };
  
  ## Additional input checks
  # Ties
  if(sum(duplicated(y))>0){warning("Rank normal transformation is not adapted for data with ties.\n")}
  # Dimension
  if(max(dim(G))>1e4){warning("Bootstrap correlation estimation will be time intensive for genotype matrices of this size.\n")}
  # Method selection
  flag.m = (method %in% c("AvgCorr","Bootstrap","Manual"));
  if(!flag.m){stop("Select 'AvgCorr', 'Bootstrap', or 'Manual' as the method for correlation estimation.")};
  
  ## Association testing
  # Calculate D-INT p-values
  P1 = suppressWarnings(DINT(y=y,G=G,X=X,S=S,M=F,c=c));
  # Calculate PI-INT p-values
  P2 = suppressWarnings(PIINT(y=y,G=G,X=X,S=S,M=F,c=c));
  # Specify correlation between z(DINT) and z(PIINT)
  if(method=="AvgCorr"){
    # Obtain correlation by averaging across loci
    R = AvgCorr(P1=P1,P2=P2);
  } else if(method=="Bootstrap") {
    # Obtain bootstrap correlation estimate
    R = BootCorr(y=y,G=G,X=X,S=S,B=B,parallel=parallel);
  } else {
    # Manually specify correlation
    if(missing(set.rho)){stop("Provide a value for rho if using method=='Manual'.")}
    R = set.rho;
  }
  # Matrix containing P1, P2, and their estimted correlation;
  Q = cbind(P1,P2,R);
  # Calculate Omnibus p-values
  POmni = apply(X=Q,MARGIN=1,FUN=Omni);
  # Output
  Out = cbind(P1,P2,POmni);
  colnames(Out) = c("DINT","PIINT","RNOmni");
  if(keep.rho){
    Out = cbind(Out,R);
    colnames(Out)[4] = "Corr";
  }
  return(Out);
}