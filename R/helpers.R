#' Smoothing Kernel for Density Estimation
#'
#' Implements Gaussian kernel smoothing.
#'
#' @param x A vector of numeric values.
#' @param mu The mean of the Gaussian density.
#' @param sigma The standard deviation of the Gaussian density.
#'
#' @return A vector of numeric values.
#' @export
#'
smoothing_kernel <- function(x, mu, sigma) {
  yseq <- stats::dnorm((x - mu) / sigma)
  return((sigma^(-length(x))) * prod(yseq))
}

#' Find the Coefficients of a Weighted Spline Function
#'
#' @param E A numeric matrix.
#' @param W A numeric matrix.
#' @param t_val A numeric matrix.
#' @param X A numeric matrix.
#' @param w The smoothing parameter.
#' @param d The intrinsic dimension.
#' @param D The dimension of the higher dimensional space.
#'
#' @return A numeric matrix.
#'
#' @noRd
solve_weighted_spline <- function(E, W, t_val, X, w, d, D) {
  M1 <- cbind(
    2 * E %*% W %*% E + 2 * w * E,
    2 * E %*% W %*% t_val,
    t_val
  )
  M2 <- cbind(
    2 * t(t_val) %*% W %*% E,
    2 * t(t_val) %*% W %*% t_val,
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M3 <- cbind(
    t(t_val),
    matrix(0, ncol = d + 1, nrow = d + 1),
    matrix(0, ncol = d + 1, nrow = d + 1)
  )
  M <- rbind(M1, M2, M3)

  b <- rbind(
    2 * E %*% W %*% X,
    2 * t(t_val) %*% W %*% X,
    matrix(0, nrow = d + 1, ncol = D)
  )
  sol <- MASS::ginv(M) %*% b
  sol
}

#' Project onto Low-Dimensional Manifold
#'
#' @param x A data point in high-dimensional space.
#' @param f An embedding function.
#' @param initial_guess Guess of the parameterization.
#'
#' @return A vector describing the data point in low-dimensional space.
#'
#' @noRd
projection_pme <- function(x, f, initial_guess) {
  est <- stats::nlm(function(t) dist_euclidean(x = x, f(t)), p = initial_guess)
  est$estimate
}
