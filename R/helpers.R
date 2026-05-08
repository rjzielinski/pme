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


#' Sample points from the Fibonacci Sphere
#'
#' @param n The number of points
#'
#' @return A matrix of points, using spherical coordinates.
#' @export
fibonacci_sphere <- function(n) {
  phi <- (1 + sqrt(5)) / 2
  points <- 1:n

  coord1 <- (points / phi) %% 1
  coord2 <- points / n

  theta <- 2 * pi * coord1
  phi <- acos(1 - 2 * coord2)

  x <- cos(theta) * sin(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(phi)

  pracma::cart2sph(cbind(x, y, z))[, -3]
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
#' @export
solve_weighted_spline_r <- function(E, W, t_val, X, w, d, D) {
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
#' @export
solve_spline_R <- function(E, t_val, X, w, d, D) {
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
      function(t) dist_euclidean(x = x, f(matrix(c(x[1], t), nrow = 1))),
      p = initial_guess[-1]
    ),
    silent = TRUE
  )
  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = initial_guess[-1],
        function(t) dist_euclidean(x = x, f(c(x[1], t))),
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(x[1], nlopt_est$solution))
    }
  } else {
    return(c(x[1], nlm_est$estimate))
  }
}

#' Project onto Low-Dimensional Manifold
#'
#' @param x A value
#' @param f A value
#' @param initial_guess A value
#'
#' @return A value
#' @export
projection_lpme_opt <- function(x, f, initial_guess) {
  x1 <- x[1]
  init_param <- initial_guess[-1]

  obj_fun <- function(t) {
    dist_euclidean(x = x, f(c(x1, t)))
  }

  nlm_est <- try(
    stats::nlm(
      obj_fun,
      p = init_param
    ),
    silent = TRUE
  )
  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = init_param,
        obj_fun,
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(x1, nlopt_est$solution))
    }
  } else {
    return(c(x1, nlm_est$estimate))
  }
}

#' Project onto Low-Dimensional Manifold
#'
#' @param x A value
#' @param f A value
#' @param initial_guess A value
#'
#' @return A value
#' @export
projection_lpme_opt2 <- function(x, f, initial_guess) {
  f_new <- function(t) {
    coefs <- f_coef_list$f(t[1])
    coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
    return_vec <- t(coef_mat[1:n_knots, ]) %*%
      etaFunc(t[-1], t_initial, 4 - d) +
      t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*%
        matrix(c(1, t[-1]), ncol = 1)
    c(t[1], return_vec)
  }

  x1 <- x[1]
  init_param <- initial_guess[-1]

  obj_fun <- function(t) {
    sum((x - f(c(x1, t)))^2)
  }

  nlm_est <- try(
    stats::nlm(
      obj_fun,
      p = init_param
    ),
    silent = TRUE
  )

  if (inherits(nlm_est, "try-error")) {
    opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-07)
    nlopt_est <- try(
      nloptr::nloptr(
        x0 = init_param,
        obj_fun,
        opts = opts
      ),
      silent = TRUE
    )
    if (inherits(nlopt_est, "try-error")) {
      return(NULL)
    } else {
      return(c(x1, nlopt_est$solution))
    }
  } else {
    return(c(x1, nlm_est$estimate))
  }
}
