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
// [[Rcpp::export]]
double norm_euclidean(arma::vec x) {
  return arma::norm(x, 2);
}

//' Compute the Euclidean Distance between Vectors
//'
//' C++ implementation of Euclidean distance calculations for use in
//' internal Rcpp functions.
//'
//' @param x Numeric vector of values.
//' @param y Numeric vector of values.
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

//' Documentation Still Needed
//'
//' This is a function that still needs to be documented properly.
//'
//' @param t Numeric vector of values.
//' @param lambda Number of dimensions.
//'
//' @return A numeric value.
// [[Rcpp::export]]
double eta_kernel(arma::vec t, int lambda) {
  double lambda_num = lambda / 1.0;
  double norm_val = norm_euclidean(t);
  double y;
  if (lambda % 2 == 0) {
    if (norm_val == 0) {
      y = 0;
    } else {
      y = pow(norm_val, lambda_num) * log(norm_val);
    }
  } else {
    y = pow(norm_val, lambda_num);
  }
  return y;
}

//' Documentation Still Needed
//'
//' This is a function that still needs to be documented properly.
//'
//' @param x Numeric matrix of values.
//' @param lambda Number of dimensions.
//'
//' @return A numeric matrix.
// [[Rcpp::export]]
arma::mat calcE(arma::mat x, int lambda) {
  int nrow = x.n_rows;
  arma::mat E(nrow, nrow);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < nrow; j++) {
      E(i, j) = eta_kernel(x.row(i).t() - x.row(j).t(), lambda);
    }
  }
  return E;
}

//' Documentation Still Needed
//'
//' This is a function that still needs to be documented properly.
//'
//' @param t A numeric vector of values.
//' @param tau A numeric matrix of values.
//' @param lambda The number of dimensions.
//'
//' @return A numeric matrix.
// [[Rcpp::export]]
arma::mat etaFunc(arma::vec t, arma::mat tau, int lambda) {
  int nrow = tau.n_rows;

  arma::mat eta(nrow, 1);

  for (int i = 0; i < nrow; i++) {
    eta(i, 0) = eta_kernel(t - tau.row(i).t(), lambda);
  }
  return eta;
}

//' Documentation Still Needed
//'
//' This is a function that still needs to be documented properly.
//'
//' @param df A numeric matrix.
//' @param x A numeric matrix.
//'
//' @return A numeric vector.
// [[Rcpp::export]]
arma::vec calc_nearest_x(arma::mat df, arma::mat x) {
  arma::vec nearest(df.n_rows);
  for (int i = 0; i < df.n_rows; i++) {
    arma::mat temp_x = x.rows(find(x.col(0) == df(i, 0)));
    arma::vec distances(temp_x.n_rows);
    for (int j = 0; j < temp_x.n_rows; j++) {
      distances(j) = dist_euclidean(df.row(i).t(), temp_x.row(j).t());
    }
    nearest(i) = distances.index_min();
  }
  return nearest;
}

//' Documentation Still Needed
//'
//' This is a function that still needs to be documented properly.
//'
//' @param df A numeric matrix.
//' @param tnew A numeric matrix.
//' @param nearest_x A numeric vector.
//'
//' @return A numeric matrix.
// [[Rcpp::export]]
arma::mat calc_init_param(arma::mat df, arma::mat tnew, arma::vec nearest_x) {
  arma::mat params(df.n_rows, tnew.n_cols);
  for (int i = 0; i < df.n_rows; i++) {
    arma::mat temp_params = tnew.rows(find(tnew.col(0) == df(i, 0)));
    params.row(i) = temp_params.row(nearest_x(i));
  }
  return params;
}
