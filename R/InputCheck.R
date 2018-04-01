# Purpose: Function to check input to association tests
# Updated: 180330

#' Input Check
#' 
#' Function to ensure the dimensions of inputs to association methods agree.
#' 
#' @param y Numeric phenotype vector.
#' @param G Snp by obs genotype matrix.
#' @param X Obs by feature covariate matrix.
#' @param S Obs by feature structure matrix.

inCheck = function(y,G,X,S){
  ## Function to check row for missingness
  aux = function(x){sum(is.na(x))==0};
  ## Check phenotype
  # Ensure continuous phenotype is supplied
  flag.y = !(is.numeric(y)&is.vector(y));
  if(flag.y){
    warning("A numeric vector is required for y.");
  }
  n.y = length(y);
  # Check for missing values
  keep = !is.na(y);
  ## Check genotype
  # Ensure entires of G are double precision
  storage.mode(G) = "double";
  # Change to matrix if vector is supplied
  if(is.vector(G)){
    G = matrix(G,nrow=1);
  }
  n.g = ncol(G);
  # Ensure numeric
  flag.g = !is.matrix(G);
  if(flag.g){
    warning("A numeric matrix is required for G.");
  }
  ## Check covariates
  # Change to data.frame if vector is supplied
  if(is.vector(X)){
    X = matrix(X,ncol=1);
  }
  n.x = nrow(X);
  # Ensure numeric
  flag.x = !is.matrix(X);
  if(flag.x){
    warning("A numeric matrix is required for X.");
  }
  # Check for missing values
  if(sum(is.na(X))>0){
    keep = keep&apply(X,MARGIN=1,FUN=aux);
    warning("Covariate matrix contains missing values.");
  }
  ## Check structure matrix
  # Change to data.frame if vector is supplied
  if(is.vector(S)){
    S = matrix(S,ncol=1);
  }
  n.s = nrow(S);
  # Ensure numeric
  flag.s = !is.matrix(S);
  if(flag.s){
    warning("A numeric matrix is required for S.");
  }
  # Check for missing values
  if(sum(is.na(S))>0){
    keep = keep&apply(X=S,MARGIN=1,FUN=aux);
    warning("Structure matrix contains missing values.");
  }
  # Dimensional consistency
  flag.d = !all.equal(n.y,n.g,n.x,n.s);
  if(flag.d){
    warning("Dimensions of inputs are inconsistent. Ensure length(y)=ncol(G)=nrow(X)=nrow(S).")
  }
  # Input failure
  fail = (flag.y|flag.x|flag.s|flag.g|flag.d);
  # Drop missing observations
  y = y[keep];
  G = G[,keep,drop=F];
  X = X[keep,,drop=F];
  S = S[keep,,drop=F];
  # Ensure y is standardized
  y = as.numeric(scale(y));
  Out = list("fail"=fail,"y"=y,"G"=G,"X"=X,"S"=S);
  return(Out);
}