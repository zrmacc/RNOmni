#' Rank-Normalize
#' 
#' Applies the rank-normal transform to a numeric vector. 
#' @param u Numeric vector.
#' @param c Offset. Defaults to (3/8), correspond to the Blom transform.
#' @return Numeric vector.
#' @importFrom stats qnorm
#' @export 

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

#' Basic Association Test
#' 
#' Regression of phenotype on genotype, covariates, and structure.
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @importFrom RcppEigen fastLmPure 
#' @importFrom stats coef model.matrix pt
#' @export

BAT = function(y,G,X,S){
  # Design matrix
  D = cbind(X,S);
  D = model.matrix(~.,data=data.frame(D));
  q = ncol(D);
  # Function to implement regression and calculate wald p-value
  getP = function(g){
    # Add genotype
    Dg = cbind(D,g);
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
#' Rank-normalizes the phenotype, then regresses the transformed
#' phenotype on genotype, covariates, and structure. 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix. 
#' @param c Offset applied during rank-normalization.
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt
#' @export
#' 
DINT = function(y,G,X,S,c=3/8){
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
#' Two-stage regression. Regresses phenotype on covariates and structure
#' to obtain residuals. Regresses transformed residuals on genotype.
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix. 
#' @param c Offset applied during rank-normalization.
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt resid
#' @export
#' 
FIINT = function(y,G,X,S,c=3/8){
  warning("This function was included for simulation purposes, and was not found to provide valid inference.");
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
#' Two-stage regression. Regresses phenotype on covariates to obtain 
#' residuals. Regresses transformed residuals on genotype and structure. 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix. 
#' @param c Offset applied during rank-normalization.
#' @importFrom RcppEigen fastLmPure
#' @importFrom stats coef model.matrix pt resid
#' @export
#' 
PIINT = function(y,G,X,S,c=3/8){
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
#' Estimate correlation using the average of z1*z2 across loci.
#' @param P1 Fist p-value.
#' @param P2 Second p-value.
#' @param eps Force correlation estimate to fall in (eps,1-eps);
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
#' @param cores Cores to use during bootstrapping if running in parallel.
#' @importFrom abind abind
#' @importFrom doMC registerDoMC
#' @importFrom stats cor qnorm
#' @import foreach

BootCorr = function(y,G,X,S,B=100,cores=1){
  # Parallelize
  if(cores>1){doMC::registerDoMC(cores=cores)};
  # Obs
  n = length(y);
  # Bind results into an array
  acomb = function(...) {abind::abind(..., along=3)};
  Rho = foreach(i=1:B,.combine="acomb",.multicombine=T,.inorder=F) %dopar% {
    # Draw sample
    Draw = sample(x=n,replace=T);
    # Bootstrap frames
    yb = y[Draw];
    Xb = X[Draw,];
    Sb = S[Draw,];
    # Calculate z-statistics
    z1 = -qnorm(DINT(y=yb,G=G,X=Xb,S=Sb));
    z2 = -qnorm(PIINT(y=yb,G=G,X=Xb,S=Sb));
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
#' Omnibus association test that estimates p-values by direct and partially 
#' indirect inverse normal transformation (INT), then synthesizes the results
#' into a single test statistic, based on whichever approach provides more
#' evidence against null hypothesis.
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.
#' @param method Method used to estimate correlation for the omnibus test,
#'   either "AvgCorr" or "Bootstrap".
#' @param c Offset applied during rank-normalization.
#' @param B Bootstrap samples for correlation estimation.
#' @param rho Logical indicating whether to return the correlation parameter
#'   estimated during omnibus calculation. Defaults to FALSE.
#' @param cores Cores to use during bootstrapping if running in parallel.
#' @export

RNOmni = function(y,G,X,S,method="AvgCorr",c=3/8,B=100,rho=F,cores=1){
  ## Input checks
  # Ties
  if(sum(duplicated(y))>0){warning("Rank normal transformation is not adapted for data with ties.")}
  # Dimension
  n = length(y);
  flag = (ncol(G)==n)&(nrow(X)==n)&(nrow(S)==n);
  if(!flag){stop("Dimensions of input data do not agree.")};
  if(max(dim(G))>1e4){warning("Bootstrap correlation estimation will be time intensive for genotype matrices of this size")}
  # Method selection
  flag = (method %in% c("AvgCorr","Bootstrap"));
  if(!flag){stop("Select 'AvgCorr' or 'Bootstrap' as the method for correlation estimation.")};
  
  ## Association testing
  # Calculate D-INT p-values
  P1 = DINT(y=y,G=G,X=X,S=S,c=c);
  # Calculate PI-INT p-values
  P2 = PIINT(y=y,G=G,X=X,S=S,c=c);
  # Branch to average correlation or bootstrap method
  if(method=="AvgCorr"){
    R = AvgCorr(P1=P1,P2=P2);
  } else {
    # Obtain bootstrap correlation estimate
    R = BootCorr(y=y,G=G,X=X,S=S,B=B,cores=cores);
  }
  # Matrix containing P1, P2, and their estimted correlation;
  Q = cbind(P1,P2,R);
  # Calculate Omnibus p-values
  POmni = apply(X=Q,MARGIN=1,FUN=Omni);
  # Output
  Out = cbind(P1,P2,POmni);
  colnames(Out) = c("DINT","PIINT","RNOmni");
  if(rho){
    Out = cbind(Out,R);
    colnames(Out)[4] = "Corr";
  }
  return(Out);
}