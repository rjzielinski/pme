#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

#include <numeric>
#include <vector>
#include <cmath>
#include <complex>

//' RcppArmadillo Hello World
//'
//' This is a simple example of creating two matrices and returning
//' the result of an operation on them.
//'
//' @export
//' @examples
//' rcpparma_hello_world
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
  arma::mat m1 = arma::eye<arma::mat>(3, 3);
  arma::mat m2 = arma::eye<arma::mat>(3, 3);

  return m1 + 3 * (m1 + m2);
}

// [[Rcpp::export]]
double norm_euclidean(arma::vec x) {
  double sq_sum = 0;
  for (int i = 0; i < x.size(); i++) {
    sq_sum += pow(x[i], 2);
  }
  return sqrt(sq_sum);
}
