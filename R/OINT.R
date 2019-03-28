# Purpose: Omnibus test
# Updated: 19/03/10

#' Omnibus p-value
#' 
#' @param p Vector of p-values
#' @return OINT p-value.
#' 
#' @importFrom stats pcauchy

OINTp = function(p){
  # Check input
  m = min(p);
  M = max(p);
  # Cases
  if((m<0)|(M>1)){stop("Cannot have p-values < 0 or > 1.")}
  if((m==0)&(M==1)){stop("Cannot have p-values of 0 and 1 simultaneously.");}
  if((m==0)&(M<1)){return(0);}
  if((m>0)&(M==1)){return(1);}
  # Convert to Cauchy
  aux = function(x){
    if(x<1e-10){
      y = 1/(x*pi);
    } else {
      y = tanpi(0.5-x);
    }
    return(y);
  }
  y = mean(sapply(p,aux));
  # Invert
  if(y>1e10){
    q = (1/y)/pi;
  } else {
    q = pcauchy(q=y,lower.tail=F);
  }
  return(q);
}

#' Omnibus-INT
#' 
#' Association test that synthesizes the \code{\link{DINT}} and
#' \code{\link{IINT}} tests. The first approach is most powerful for traits that
#' could have arisen from a rank-preserving transformation of a latent normal
#' trait. The second approach is most powerful for traits that are linear in
#' covariates, yet have skewed or kurtotic residual distributions. During the
#' omnibus test, the direct and indirect tests are separately applied then
#' 
#' @importFrom stats qnorm
#' @importFrom plyr aaply
#' @export
#' 
#' @param y Numeric phenotype vector.
#' @param G Obs by locus genotype matrix.
#' @param X Model matrix of covariates and structure adjustments. Should include
#'   an intercept. Omit to perform marginal tests of association.
#' @param k Offset applied during rank-normalization. See 
#'   \code{\link{rankNorm}}.
#' @param simple Return the OINT p-values only? 
#' @param parallel Logical indicating whether to run in parallel. Must register 
#'   parallel backend first.
#' @return A numeric matrix of p-values, three for each column of \code{G}.
#'   
#' @seealso Basic association test \code{\link{BAT}}, direct INT \code{\link{DINT}}, 
#' indirect INT \code{\link{IINT}}.
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
#' # Omnibus
#' p = OINT(y=y,G=G,X=X,simple=T);
#' }

OINT = function(y,G,X=NULL,k=3/8,simple=FALSE,parallel=FALSE){
  ## Check inputs
  # Input check 
  n = length(y);
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(G)){stop("A numeric matrix is expected for G.")};
  ng = ncol(G);
  if(is.null(X)){X=array(1,dim=c(n,1))};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Please exclude observations missing phenotype or covariate information.")}

  ## Association testing
  # Calculate D-INT p-values
  p1 = DINT(y=y,G=G,X=X,k=k,parallel=parallel,simple=T);
  # Calculate PI-INT p-values
  p2 = IINT(y=y,G=G,X=X,k=k,parallel=parallel,simple=T);
  
  # P Matrix
  P = cbind(p1,p2);
  
  # Omnibus p-values
  p3 = aaply(.data=P,.margins=1,.fun=OINTp);
  
  ## Output
  
  # Locus names
  gnames = colnames(G);
  if(is.null(gnames)){
    gnames = seq(1:ng);
  }
  
  # Format
  if(simple){
    Out = p3;
    names(Out) = gnames;
  } else {
    Out = cbind("DINT-p"=p1,"IINT-p"=p2,"OINT-p"=p3);
    rownames(Out) = gnames;
  }
  
  # Return
  return(Out);
}