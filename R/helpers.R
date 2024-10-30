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


#' Smoothing Kernel for Density Estimation
#' 
#' Implements Gaussian kernel smoothing on log-scale
#' 
#' @param x A vector of numeric values.
#' @param mu The mean of the Gaussian density
#' @param sigma The standard deviation of the Gaussian density
#' 
#' @return A numeric value
#' @export
#' 
log_smoothing_kernel_r <- function(x, mu, sigma) {
  yseq <- stats::dnorm((x - mu) / sigma)
  output <- (-length(x) * log(sigma)) + sum(log(yseq))
  return(output)
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

#' Solve for Smoothing Spline Coefficients
#'
#' @param E A matrix of values.
#' @param t_val A numeric matrix of input values.
#' @param X A numeric matrix of outputs.
#' @param w The smoothing parameter.
#' @param d The dimension of the input.
#' @param D The dimension of the output.
#'
#' @return A numeric matrix of coefficients.
#'
#' @noRd
solve_spline <- function(E, t_val, X, w, d, D) {
  M1 <- cbind(E + (w * diag(rep(1, nrow(t_val)))), t_val)
  M2 <- cbind(t(t_val), matrix(0, ncol = d + 1, nrow = d + 1))
  M <- rbind(M1, M2)
  b <- rbind(X, matrix(0, nrow = d + 1, ncol = D))
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
#' @export
projection_pme <- function(x, f, initial_guess) {
  nlm_est <- try(
    stats::nlm(
      function(t) dist_euclidean(x = x, f(t)),
      p = initial_guess
    ),
    silent = TRUE
  )

  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-10)
     nlopt_est <- try(
       nloptr::nloptr(
        x0 <- initial_guess,
        function(t) dist_euclidean(x = x, f(t)),
        opts = opts
      ),
      silent = TRUE
     )
     if (inherits(nlopt_est, "try-error")) {
       return(NULL)
     } else {
       return(nlopt_est$solution)
     }
  } else {
    return(nlm_est$estimate)
  }
}

#' Project onto Low-Dimensional Manifold
#'
#' @param x A value
#' @param f A value
#' @param initial_guess A value
#' @param n_knots A value
#' @param d_new A value
#' @param gamma A value
#'
#' @return A value
#' @export
projection_lpme <- function(x, f, initial_guess, n_knots, d_new, gamma) {
  nlm_est <- try(
    stats::nlm(
      function(t) dist_euclidean(x = x, f(matrix(c(initial_guess[1], t), nrow = 1))),
      p = initial_guess[-1]
    ),
    silent = TRUE
  )
  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = initial_guess[-1],
        function(t) dist_euclidean(x = x, f(c(initial_guess[1], t))),
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(initial_guess[1], nlopt_est$solution))
    }
  } else {
    return(c(initial_guess[1], nlm_est$estimate))
  }
}
