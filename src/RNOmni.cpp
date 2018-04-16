// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for calculating score tests in RNOmni models
// Updated: 180330
#include <RcppEigen.h>

//' Normal Model
//' 
//' @param y Outcome.
//' @param Z Model matrix.
// [[Rcpp::export]]

SEXP fitNorm(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Z){
  // Observations
  const int n = y.size();
  // Estimated parameters
  const int p = Z.cols();
  // Gram matrix
  const Eigen::MatrixXd ZtZ = Z.transpose()*Z;
  // Estimate beta
  const Eigen::VectorXd b = (ZtZ).llt().solve(Z.transpose()*y);
  // Calculate residuals
  const Eigen::VectorXd eT = (y-Z*b);
  // Scale
  const double qf = (eT.transpose()*eT);
  const double tau = qf/(n-p);
  // Information
  const Eigen::MatrixXd Ibb = ZtZ/tau;
  return Rcpp::List::create(Rcpp::Named("Beta")=b,Rcpp::Named("Tau")=tau,Rcpp::Named("Ibb")=Ibb,Rcpp::Named("eT")=eT);
}

//' Correlation
//' 
//' Calculates the correlation between two vectors.
//' 
//' @param a First vector.
//' @param b Second vector.
// [[Rcpp::export]]

SEXP vecCor(const Eigen::Map<Eigen::VectorXd> a,const Eigen::Map<Eigen::VectorXd> b){
  const int n = a.size();
  // One vector
  const Eigen::VectorXd j = Eigen::VectorXd::Constant(n,1);
  // Centering
  const double at = j.transpose()*a;
  const double am = at/n;
  const Eigen::VectorXd ac = a-am*j;
  const double bt = j.transpose()*b;
  const double bm = bt/n;
  const Eigen::VectorXd bc = b-bm*j;
  // Correlation
  const double acbc = ac.transpose()*bc;
  const double acac = ac.transpose()*ac;
  const double bcbc = bc.transpose()*bc;
  const double r = acbc/sqrt(acac*bcbc);
  return Rcpp::wrap(r);
} 

