#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>
#include <boost/math/distributions/normal.hpp>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]


using namespace Rcpp;
using namespace std;
using boost::math::normal;

#include <numeric>
#include <vector>
#include <cmath>
#include <complex>
#include <random>

//' Compute the Euclidean Norm of a Vector
//'
//' C++ implementation of the Euclidean norm calculations for a vector. This
//' is meant for internal use for the Rcpp implementation.
//'
//' @param x A numeric vector of values.
//'
//' @return A numeric value.
//' @export
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

//' Sum in Logspace
//'
//' Given log(x) and log(y), return log(x+y)
//'
//' @param lx A double value in logspace
//' @param ly A double value in logspace
//' @return A numeric value
//' @export
// [[Rcpp::export]]
double logspace_sum(double lx, double ly) {
  double log_max = max(lx, ly);
  double log_sum = log_max + log(1 + exp(-abs(lx - ly)));
  return(log_sum);
}

//' Difference in Logspace
//'
//' Given log(x) and log(y), return log(x-y)
//'
//' @param lx A double value in logspace
//' @param ly A double value in logspace
//' @return A numeric value
//' @export
// [[Rcpp::export]]
double logspace_diff(double lx, double ly) {
  // double log_diff = lx + log(1 - exp(ly - lx));
  double log_max = max(lx, ly);
  double log_min = min(lx, ly);
  double log_diff = log_max + log(1 - exp(log_min - log_max));
  return(log_diff);
}

//' Sum over a Vector in Logspace
//'
//' Recursively call logspace_sum() over n > 2 values
//'
//' @param x A numeric vector of values in logspace
//' @return A numeric value
//' @export
//'
// [[Rcpp::export]]
double logspace_sum_vec(arma::vec x) {
  int n = x.n_elem;
  double log_sum;
  if (n >= 2) {
    log_sum = logspace_sum(x(0), x(1));
    if (n > 2) {
      for (int i = 2; i < n; i++) {
        log_sum = logspace_sum(log_sum, x(i));
      }
    }
  } else {
    throw std::range_error("x must have length >= 2.");
  }
  return(log_sum);
}

//' Smoothing Kernel for Density Estimation
//'
//' Implements Gaussian kernel smoothing on log-scale
//'
//' @param x A numeric vector of values.
//' @param mu The mean of the Gaussian density
//' @param sigma The standard deviation of the Gaussian density
//'
//' @return A numeric value.
//' @export
// [[Rcpp::export]]
double log_smoothing_kernel(arma::vec x, arma::vec mu, double sigma) {
  //NumericVector x_mod = as<NumericVector>(wrap(x));
  //NumericVector mu_mod = as<NumericVector>(wrap(mu));
  //NumericVector z_probs = dnorm((x_mod - mu_mod) / sigma);
  //int n = x_mod.length();
  int n = x.n_elem;
  NumericVector z = as<NumericVector>(wrap((x - mu) / sigma));
  NumericVector z_probs = dnorm(z);
  double output = (-n * log(sigma)) + sum(log(z_probs));
  // Avoid outputting infinite values for extremely small probabilities
  if (isinf(output)) {
    output = -DBL_MAX;
  }
  return(output);
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
//' @export
//'
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
//' @export
//'
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

//' Calculate "A" Matrix
//'
//' This is a function that still needs to be documented properly.
//'
//' @param x_obs A numeric matrix containing unreduced data.
//' @param mu A numeric matrix of component centers.
//' @param sigma A numeric value denoting the bandwidth of density estimation.
//'
//' @return A numeric matrix.
// [[Rcpp::export]]
arma::mat calc_A(arma::mat x_obs, arma::mat mu, double sigma) {
  int n = x_obs.n_rows;
  int N = mu.n_rows;
  arma::mat A(n, N);
  x_obs = x_obs.t();
  mu = mu.t();
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < n; i++) {
      A(i, j) = log_smoothing_kernel(x_obs.col(i), mu.col(j), sigma);
    }
  }
  return A;
}

//' Calculate Weights of Mixture Components
//'
//' @param x_obs A numeric matrix containing the unreduced data.
//' @param mu A numeric matrix of component centers.
//' @param sigma A numeric value denoting the bandwidth of the density estimation.
//' @param epsilon A numeric value denoting the tolerance of the Euclidean distance between weights.
//' @param max_iter An integer denoting the maximum number of iterations.
//'
//' @return A numeric vector of weights.
// [[Rcpp::export]]
arma::vec calc_weights_cpp(arma::mat x_obs, arma::mat mu, double sigma, double epsilon, int max_iter) {
  int n = x_obs.n_rows;
  int D = x_obs.n_cols;
  int N = mu.n_rows;
  arma::mat A = calc_A(x_obs, mu, sigma);
  return(A.row(0));
}

//calc_weights <- function(x_obs, mu, sigma, epsilon = 0.001, max_iter = 1000) {
//  # Initialize function parameters
//  n <- nrow(x_obs)
//  D <- ncol(x_obs)
//  N <- nrow(mu)
//  A <- calc_A(x_obs, mu, sigma)
//  W <- t(t(A) + log(theta_old))
//  theta_old <- rep(1 / N, N)
//  abs_diff <- 10 * epsilon
//  count <- 0
//  lambda_hat_old <- c(n, rep(-1, D))


//  while ((abs_diff > epsilon) & (count <= max_iter)) {
//    # W <- t(t(A) * theta_old)
//    W <- t(t(A) + log(theta_old))
//    w <- Rfast::colsums(W / Rfast::rowsums(W))
//    lambda_hat <- stats::nlm(
//      f = f_lambda,
//      p = lambda_hat_old,
//      x = x_obs,
//      mu = mu,
//      w = w,
//      iterlim = 1000
//    )$estimate
//    theta_new <- w / Rfast::rowsums(t(t(cbind(rep(1, N), mu)) * lambda_hat))

//    abs_diff <- dist_euclidean(theta_new, theta_old)
//    if (is.na(abs_diff)) {
//      return(bound_theta(theta_old))
//    } else {
//      theta_old <- bound_theta(theta_new)
//      count <- count + 1
//      lambda_hat_old <- lambda_hat
//    }
//  }
//  bound_theta(theta_new)
//}

//' Find the Coefficients of a Weighted Spline Function
//'
//' @param E A numeric matrix.
//' @param W A numeric matrix.
//' @param t_val A numeric matrix.
//' @param X A numeric matrix.
//' @param w The smoothing parameter.
//' @param d The intrinsic dimension.
//' @param D The dimension of the higher dimensional space.
//'
//' @return A numeric matrix.
arma::mat solve_weighted_spline_reg(arma::mat E, arma::mat W, arma::mat t_val, arma::mat X, double w, int d, int D, double jitter = 1e-8) {
  arma::mat M1 = join_rows(2 * E * W * E + (2 * w * E), 2 * E * W * t_val);
  M1 = join_rows(M1, t_val);
  arma::mat M2 = join_rows(2 * t_val.t() * W * E, 2 * t_val.t() * W * t_val);
  arma::mat zero_mat = arma::zeros(d + 1, d + 1);
  M2 = join_rows(M2, zero_mat);
  arma::mat M3 = join_rows(t_val.t(), zero_mat);
  M3 = join_rows(M3, zero_mat);
  arma::mat M = join_cols(M1, M2);
  M = join_cols(M, M3);
  arma::mat b = join_cols(2 * E * W * X, 2 * t_val.t() * W * X);
  arma::mat zero_mat2 = arma::zeros(d + 1, D);
  b = join_cols(b, zero_mat2);
  // arma::mat sol = arma::pinv(M) * b;
  // for computational efficiency, use arma::solve() instead of Moore-Penrose pseudoinverse
  // M is often singular, so approximate the solution by adding small jitter
  M.diag() += jitter;
  arma::mat sol = arma::solve(M, b);
  return sol;
}

//' Find the Coefficients of a Weighted Spline Function
//'
//' @param E A numeric matrix.
//' @param W A numeric matrix.
//' @param t_val A numeric matrix.
//' @param X A numeric matrix.
//' @param w The smoothing parameter.
//' @param d The intrinsic dimension.
//' @param D The dimension of the higher dimensional space.
//'
//' @return A numeric matrix.
arma::mat solve_weighted_spline_pinv(arma::mat E, arma::mat W, arma::mat t_val, arma::mat X, double w, int d, int D, double jitter = 1e-8) {
  arma::mat M1 = join_rows(2 * E * W * E + (2 * w * E), 2 * E * W * t_val);
  M1 = join_rows(M1, t_val);
  arma::mat M2 = join_rows(2 * t_val.t() * W * E, 2 * t_val.t() * W * t_val);
  arma::mat zero_mat = arma::zeros(d + 1, d + 1);
  M2 = join_rows(M2, zero_mat);
  arma::mat M3 = join_rows(t_val.t(), zero_mat);
  M3 = join_rows(M3, zero_mat);
  arma::mat M = join_cols(M1, M2);
  M = join_cols(M, M3);
  arma::mat b = join_cols(2 * E * W * X, 2 * t_val.t() * W * X);
  arma::mat zero_mat2 = arma::zeros(d + 1, D);
  b = join_cols(b, zero_mat2);
  arma::mat sol = arma::pinv(M) * b;
  // for computational efficiency, use arma::solve() instead of Moore-Penrose pseudoinverse
  // M is often singular, so approximate the solution by adding small jitter
  // M.diag() += jitter;
  // arma::mat sol = arma::solve(M, b);
  return sol;
}

//' Find the Coefficients of a Weighted Spline Function
//'
//' @param E A numeric matrix.
//' @param W A numeric matrix.
//' @param t_val A numeric matrix.
//' @param X A numeric matrix.
//' @param w The smoothing parameter.
//' @param d The intrinsic dimension.
//' @param D The dimension of the higher dimensional space.
//'
//' @return A numeric matrix.
//' @export
// [[Rcpp::export]]
arma::mat solve_weighted_spline(arma::mat E, arma::mat W, arma::mat t_val, arma::mat X, double w, int d, int D) {
  arma::mat M1 = E;
  M1.diag() += w / W.diag();
  M1 = join_rows(M1, t_val);
  arma::mat M2 = join_rows(t_val.t(), arma::zeros(d + 1, d + 1));
  arma::mat M = join_cols(M1, M2);
  arma::mat b = join_cols(X, arma::zeros(d + 1, D));
  // arma::mat sol = arma::pinv(M) * b;
  // for computational efficiency, use arma::solve() instead of Moore-Penrose pseudoinverse
  // M is often singular, so approximate the solution by adding small jitter
  // M.diag() += jitter;
  arma::mat sol = arma::solve(M, b);
  
  return sol;
}

//' Find the Coefficients of a Weighted Spline Function
//'
//' @param E A numeric matrix.
//' @param W A numeric matrix.
//' @param t_val A numeric matrix.
//' @param X A numeric matrix.
//' @param w The smoothing parameter.
//' @param d The intrinsic dimension.
//' @param D The dimension of the higher dimensional space.
//'
//' @return A numeric matrix.
arma::mat solve_weighted_spline_reinsch_pinv(arma::mat E, arma::mat W, arma::mat t_val, arma::mat X, double w, int d, int D, double jitter = 1e-8) {
  arma::mat M1 = E;
  M1.diag() += w / W.diag();
  M1 = join_rows(M1, t_val);
  arma::mat M2 = join_rows(t_val.t(), arma::zeros(d + 1, d + 1));
  arma::mat M = join_cols(M1, M2);
  arma::mat b = join_cols(X, arma::zeros(d + 1, D));
  arma::mat sol = arma::pinv(M) * b;
  // for computational efficiency, use arma::solve() instead of Moore-Penrose pseudoinverse
  // M is often singular, so approximate the solution by adding small jitter
  // M.diag() += jitter;
  // arma::mat sol = arma::solve(M, b);
  return sol;
}

