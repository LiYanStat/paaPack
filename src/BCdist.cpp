#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double BCdistcpp(arma::vec X, arma::vec Y) {
  // function to compute the BC distance for pair of amalgamation
  int n = X.n_elem;
  arma::vec Z = X + Y;
  double out = 0;
  for(int i = 0;i < n-1;i++) {
    for(int j = i+1;j < n;j++){
      double temp = as_scalar(min(X.row(i), X.row(j)) + min(Y.row(i), Y.row(j)) - min(Z.row(i), Z.row(j)));
      out += temp * temp;
    }
  }
  return out;
}