#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

#include <numeric>
#include <vector>
#include <cmath>
#include <complex>

//' Compute the Euclidean Norm of a Vector
//'
//' C++ implementation of the Euclidean norm calculations for a vector. This
//' is meant for internal use for the Rcpp implementation.
//'
//' @param x A numeric vector of values.
//'
//' @return A numeric value.
//' @export
//'
//' @examples
//' x <- 1:10
//' norm_euclidean(x)
// [[Rcpp::export]]
double norm_euclidean(arma::vec x) {
  return arma::norm(x, 2);
}

//' Compute the Euclidean Distance between Vectors
//'
//' C++ implementation of Euclidean distance calculations for use in
//' internal Rcpp functions.
//'
//' @param x, y Numeric vectors of values.
//'
//' @return A numeric value.
//' @export
//'
//' @examples
//' x <- 1:10
//' y <- 91:100
//' dist_euclidean(x, y)
// [[Rcpp::export]]
double dist_euclidean(arma::vec x, arma::vec y) {
  return norm_euclidean(x - y);
}
