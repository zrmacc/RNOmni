# Purpose: Indirect INT-based methods
# Updated: 180501

#' Indirect-INT, Method A
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates and structure adjustments to obtain residuals. In the second stage, 
#' INT-transformed residuals are regressed on genotype and population structure.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param S Model matrix of structure adjustments.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first.
#' @param check Logical indicating whether to check the input.
#' @return A numeric matrix of score statistics, one for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{calcP=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples 
#' # IINTa versus a normal phenotype
#' p = RNOmni::IINTa(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

IINTa = function(y,G,X,S,calcP=T,k=3/8,parallel=F,check=T){
  if(check){
    # Check inputs
    Input = inCheck(y,G,X,S);
    if(Input$fail){stop("Input check failed.")};
    y = Input$y;
    G = Input$G;
    X = Input$X;
    S = Input$S;
  }
  # Loci
  n.g = ncol(G);
  # Degrees of freedom
  n = length(y);
  p = ncol(X);
  q = ncol(S);
  df2 = n-q;
  # Stage 1 model
  M1 = fitNorm(y=y,Z=cbind(X,S));
  e = rankNormal(u=M1$Resid);
  # Stage 2 model
  M2 = fitNorm(y=e,Z=S);
  # Extract components
  d = M2$Resid;
  tau = M2$Tau; 
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    S.obs = S[keep,,drop=F];
    d.obs = d[keep];
    # Information components
    I11 = sum(g.obs^2);
    I12 = fastIP(A=g.obs,B=S.obs);
    I22 = fastIP(S.obs,S.obs);
    # Score statistic
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    a = as.numeric(fastIP(A=g.obs,B=d.obs));
    Ts = a^2/(tau*V);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(U,nrow=n.g);
  colnames(Out) = "Score";
  # Calculate p values
  if(calcP){
    #P = pf(q=U,df1=1,df2=df2,lower.tail=F);
    P = pchisq(q=U,df=1,lower.tail=F)
    Out = cbind(Out,P);
  }
  return(Out);
}

#' Indirect-INT, Method B
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates to obtain residuals. In the second stage, INT-transformed residuals are regressed on genotype 
#' and structure adjustments, having projected away dependence on X.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq var
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param S Model matrix of structure adjustments.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @param check Logical indicating whether to check the input.
#' @return A numeric matrix of Wald statistics, one for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{calcP=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # IINTb against normal phenotype 
#' p = RNOmni::IINTb(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

IINTb = function(y,G,X,S,calcP=T,k=3/8,parallel=F,check=T){
  if(check){
    # Check inputs
    Input = inCheck(y,G,X,S);
    if(Input$fail){stop("Input check failed.")};
    y = Input$y;
    G = Input$G;
    X = Input$X;
    S = Input$S;
  }
  # Loci
  n.g = ncol(G);
  # Degrees of freedom
  n = length(y);
  p = ncol(X);
  q = ncol(S);
  df2 = n-q;
  # Stage 1 model
  M1 = fitNorm(y=y,Z=X);
  e = rankNormal(u=M1$Resid);
  E = Resid(X=X,Y=S);
  # Stage 2 model
  M2 = fitNorm(y=e,Z=E);
  # Extract components
  d = M2$Resid;
  tau = M2$Tau; 
  # Function to calculate score statistics
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    E.obs = E[keep,,drop=F];
    d.obs = d[keep];
    # Information components
    I11 = sum(g.obs^2);
    I12 = fastIP(A=g.obs,B=E.obs);
    I22 = fastIP(A=E.obs,B=E.obs);
    # Score statistic
    V = as.numeric(SchurC(I11=I11,I22=I22,I12=I12));
    a = as.numeric(fastIP(A=g.obs,B=d.obs));
    Ts = a^2/(tau*V);
    return(Ts);
  }
  # Calculate score statistics
  U = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(U,nrow=n.g);
  colnames(Out) = "Score";
  # Calculate p values
  if(calcP){
    #P = pf(q=U,df1=1,df2=df2,lower.tail=F);
    P = pchisq(q=U,df=1,lower.tail=F)
    Out = cbind(Out,P);
  }
  return(Out);
};

#' Indirect-INT, Method C
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates and structure adjustments to obtain residuals. In the second stage,
#' INT-transformed residuals are regressed on genotype only.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq var
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param S Model matrix of structure adjustments.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @param check Logical indicating whether to check the input.
#' @return A numeric matrix of Wald statistics, one for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{calcP=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # IINTc against normal phenotype 
#' p = RNOmni::IINTc(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

IINTc = function(y,G,X,S,calcP=T,k=3/8,parallel=F,check=T){
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
  # Degrees of freedom
  df2 = length(y)-ncol(Z);
  # Fit null model
  M0 = fitNorm(y=y,Z=Z);
  # Transformed Residuals
  e = rankNormal(M0$Resid);
  # Calculate F statistic
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    e.obs = e[keep];
    # Information component
    g2 = sum(g.obs^2);
    # Wald statistic
    r2 = as.numeric(fastIP(A=g.obs,B=e.obs))^2;
    Tw = r2/(g2);
    return(Tw);
  }
  # Wald statistics
  W = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(W,nrow=n.g);
  colnames(Out) = "Wald";
  # Calculate p values
  if(calcP){
    P = pchisq(q=W,df=1,lower.tail=F);
    Out = cbind(Out,P);
  }
  return(Out);
};

#' Indirect-INT, Method D
#' 
#' Two-stage regression procedure. In the first stage, phenotype is regressed on
#' covariates and structure adjustments to obtain residuals. In the second stage,
#' INT-transformed residuals are regressed on genotype, having projected away
#' dependence on covariates and structure adjustments.
#' 
#' @importFrom plyr aaply
#' @importFrom stats pchisq var
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by snp genotype matrix.
#' @param X Model matrix of covariates.
#' @param S Model matrix of structure adjustments.
#' @param calcP Logical indicating that p values should be calculated.
#' @param k Offset applied during rank-normalization. See
#'   \code{\link{rankNormal}}.
#' @param parallel Logical indicating whether to run in parallel. Must register
#'   parallel backend first. 
#' @param check Logical indicating whether to check the input.
#' @return A numeric matrix of Wald statistics, one for each locus in \code{G},
#'   assessing the null hypothesis that genotype is unrelated to the outcome. If
#'   \code{calcP=T}, a p-value is additionally calculated for each locus.
#'   
#' @examples
#' # IINTc against normal phenotype 
#' p = RNOmni::IINTd(y=RNOmni::Y[,1],G=RNOmni::G[,1:10],X=RNOmni::X,S=RNOmni::S);

IINTd = function(y,G,X,S,calcP=T,k=3/8,parallel=F,check=T){
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
  # Degrees of freedom
  df2 = length(y)-ncol(Z);
  # Fit null model
  M0 = fitNorm(y=y,Z=Z);
  # Transformed Residuals
  e = rankNormal(M0$Resid);
  # Calculate F statistic
  aux = function(g){
    # Adjust for missingness
    keep = !is.na(g);
    g.obs = g[keep];
    e.obs = e[keep];
    # Projected genotype
    Z.obs = Z[keep,,drop=F];
    g.obs = Resid(X=Z.obs,Y=g.obs);
    # Information component
    g2 = sum(g.obs^2);
    # Wald statistic
    r2 = as.numeric(fastIP(A=g.obs,B=e.obs))^2;
    Tw = r2/(g2);
    return(Tw);
  }
  # Wald statistics
  W = aaply(.data=G,.margins=2,.fun=aux,.parallel=parallel);
  # Output frame
  Out = matrix(W,nrow=n.g);
  colnames(Out) = "Wald";
  # Calculate p values
  if(calcP){
    P = pchisq(q=W,df=1,lower.tail=F);
    Out = cbind(Out,P);
  }
  return(Out);
};