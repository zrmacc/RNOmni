// [[Rcpp::depends(RcppEigen)]]
// Purpose: Functions for calculating score tests in RNOmni models
// Updated: 180330
#include <RcppEigen.h>

//' Error Projection
//' 
//' @param X Model matrix
//' @export 
// [[Rcpp::export]]

SEXP errProj(const Eigen::Map<Eigen::MatrixXd> X){
  // Observations
  const int n = X.rows();
  // Identity
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);   
  // Projection
  const Eigen::MatrixXd P = X*(X.transpose()*X).llt().solve(X.transpose());
  const Eigen::MatrixXd Q = (I-P);
  // Output
  return Rcpp::wrap(Q);
}

//' Scale Parameter
//' 
//' @param y Outcome
//' @param Q Error projection
//' @param df Degrees of freedom
//' @export 
// [[Rcpp::export]]

SEXP scaleParam(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Q,
                const int df){
  // Quadratic form
  const double q = y.transpose()*Q*y;
  // Scale
  const double s = q/df;
  // Output
  return Rcpp::wrap(s);
}

//' Score Statistic
//' 
//' @param y Outcome
//' @param Q Error projection
//' @param g Genotype vector
//' @param s2 Scale parameter
//' @export 
// [[Rcpp::export]]

SEXP scoreStat(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Q,
               const Eigen::Map<Eigen::VectorXd> g, const double s2){
  // Variance
  const double qf = g.transpose()*Q*g;
  const double v = s2*qf;
  // Score
  const double sc = g.transpose()*Q*y;
  // Score statistic
  const double t = (sc*sc)/v;
  // Output
  return Rcpp::wrap(t);
}

//' Residuals
//' 
//' @param y Outcome
//' @param Q Error projection
//' @export 
// [[Rcpp::export]] 

SEXP eps(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::MatrixXd> Q){
  // Residual
  const Eigen::VectorXd e = Q*y;
  // Output
  return Rcpp::wrap(e);
}

//' Dot Product
//' 
//' @param a First vector
//' @param b Second vector
//' @export 
// [[Rcpp::export]] 

SEXP dotP(const Eigen::Map<Eigen::VectorXd> a, const Eigen::Map<Eigen::VectorXd> b){
  // Residual
  const Eigen::VectorXd Out = a.transpose()*b;
  // Output
  return Rcpp::wrap(Out);
}