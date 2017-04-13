#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

#include <RcppEigen.h>

using namespace Rcpp;
// using Eigen::MatrixXd;
// using namespace RcppEigen;

// [[Rcpp::export]]

SEXP eigen_mat_mult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}
