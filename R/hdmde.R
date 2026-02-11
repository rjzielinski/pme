#' Create New HDMDE Object
#'
#' @param mu A numeric matrix containing the density component centers.
#' @param sigma A numeric value describing the cluster variance.
#' @param theta A numeric vector of weight values.
#' @param km An object output from the `kmeans()` function.
#'
#' @return An object of type hdmde.
#'
#' @noRd
new_hdmde <- function(mu, sigma, theta, km) {
  hdmde_list <- list(
    mu = mu,
    sigma = sigma,
    theta_hat = theta,
    km = km
  )
  vctrs::new_vctr(hdmde_list, class = "hdmde")
}

#' Check Whether an Object is of Type `hdmde`
#'
#' @param x An object to be checked.
#'
#' @return A logical value.
#'
#' @noRd
is_hdmde <- function(x) {
  inherits(x, "hdmde")
}

#' High-Dimensional Mixture Density Estimation
#'
#' Data reduction approach to estimate a high-dimensional dataset using a
#' mixture of distributions. The centers of these distributions are identified
#' using k-means clustering.
#'
#' @param x_obs A numeric matrix containing the data to be reduced.
#' @param N0 An integer specifying the lower bound for the number of density components.
#' @param alpha A numeric value between 0 and 1 specifying the confidence level.
#' @param max_comp An integer specifying the upper bound for the number of density components.
#'
#' @return An object of type hdmde
#' @export
hdmde <- function(x_obs, N0, alpha, max_comp) {
  # Initialization ----------------------------
  zalpha <- stats::qnorm(1 - (alpha / 2))
  n <- nrow(x_obs)
  D <- ncol(x_obs)
  N <- N0

  component_estimates <- compute_estimates(x_obs, N)

  p_old <- purrr::map(
    1:n,
    ~ f_test(
      x_obs[.x, ],
      component_estimates$mu,
      component_estimates$sigma,
      component_estimates$theta_hat
    )
  ) %>%
    unlist()

  test_rejection <- TRUE

  while ((test_rejection == TRUE) & (N < min(n - 1, max_comp))) {
    N <- N + 1
    components_new <- compute_estimates(x_obs, N)
    p_new <- purrr::map(
      1:n,
      ~ f_test(
        x_obs[.x, ],
        components_new$mu,
        components_new$sigma,
        components_new$theta_hat
      )
    ) %>%
      unlist()

    difference <- p_new - p_old
    mean_diff <- mean(difference)
    sigma_hat_sq <- mean((difference - mean_diff)^2)
    Z_I_N <- sqrt(n) * mean(difference) / sqrt(sigma_hat_sq)

    if (!is.na(Z_I_N) & (abs(Z_I_N) <= zalpha)) {
      test_rejection <- FALSE
    }
    p_old <- p_new
  }

  new_hdmde(
    components_new$mu,
    components_new$sigma,
    components_new$theta_hat,
    components_new$km
  )
}

# HDMDE HELPER FUNCTIONS -------------------------------------------------------

#' Compute Mixture Component Estimates
#'
#' For a given value of N, estimate the mixture component centers, variance, and weights.
#'
#' @param x A numeric matrix.
#' @param n The number of mixture components to be estimated.
#'
#' @return A list containing the mixture centers, variance, and weights.
#'
#' @noRd
compute_estimates <- function(x, n) {
  # increasing iter.max helps to avoid poor fit in high dimensional situations.
  # km <- stats::kmeans(x, n, iter.max = 10000, nstart = 100, algorithm = "Lloyd")
  km <- stats::kmeans(x, n, iter.max = 10000, nstart = 1000)
  mu <- km$centers
  sigma_est <- estimate_sigma(x, km)
  theta_hat <- calc_weights(x, mu, sigma_est)

  list(
    mu = mu,
    sigma = sigma_est,
    theta_hat = theta_hat,
    km = km
  )
}

#' Estimate Mixture Component Variance
#'
#' @param x A numeric matrix containing the dataset.
#' @param km A k-means object including the estimated mixture components.
#'
#' @return A numeric estimate of the variance of the mixture components.
#'
#' @noRd
estimate_sigma <- function(x, km) {
  sigma_vec <- km$withinss / km$size
  sqrt(mean(sigma_vec) / ncol(km$centers))
}

#' Calculate Weights of Mixture Components
#'
#' @param x_obs A numeric matrix containing the unreduced data.
#' @param mu A numeric matrix of component centers.
#' @param sigma A numeric value denoting the bandwidth of the density estimation.
#' @param epsilon A numeric value denoting the tolerance of the Euclidean distance between weights.
#' @param max_iter The maximum number of iterations.
#'
#' @return A numeric vector of weights.
#'
#' @noRd
calc_weights <- function(x_obs, mu, sigma, epsilon = 0.001, max_iter = 1000) {
  # Initialize function parameters
  n <- nrow(x_obs)
  D <- ncol(x_obs)
  N <- nrow(mu)
  A <- exp(calc_A(x_obs, mu, sigma))

  theta_old <- rep(1 / N, N)
  abs_diff <- 10 * epsilon
  count <- 0
  lambda_hat_old <- c(n, rep(-1, D))

  while ((abs_diff > epsilon) & (count <= max_iter)) {
    W <- t(t(A) * theta_old)
    W_adj <- W / Rfast::rowsums(W)
    w <- Rfast::colsums(W_adj)

    lambda_hat <- stats::nlm(
      f = f_lambda,
      p = lambda_hat_old,
      x = x_obs,
      mu = mu,
      w = w,
      iterlim = 10000
    )$estimate
    theta_new <- w / Rfast::rowsums(t(t(cbind(rep(1, N), mu)) * lambda_hat))

    abs_diff <- max(abs(theta_new - theta_old))
    if (is.na(abs_diff) | sum(bound_theta(theta_new) == 0)) {
      return(bound_theta(theta_old))
    } else {
      theta_old <- bound_theta(theta_new)
      count <- count + 1
      lambda_hat_old <- lambda_hat
    }
  }
  bound_theta(theta_new)
}

#' HDMDE F-Test
#'
#' Hypothesis test to determine whether additional density components are necessary.
#'
#' @param x
#'
#' @return A numeric value
#'
#' @noRd
f_test <- function(x, mu, sigma, theta_hat) {
  comp_vec <- purrr::map(
    seq_len(nrow(mu)),
    ~ smoothing_kernel(x, mu[.x, ], sigma)
  ) %>%
    unlist()
  sum(theta_hat * comp_vec)
}


# CALC_WEIGHTS HELPER FUNCTIONS -----------------------------------------------

#' Calculate `A` Matrix
#'
#' @param x_obs A numeric matrix containing the unreduced dataset.
#' @param mu A numeric matrix of the component centers.
#' @param sigma A numeric value indicating the kernel bandwidth.
#'
#' @return A numeric matrix.
#'
#' @noRd
calc_A_r <- function(x_obs, mu, sigma) {
  n <- nrow(x_obs)
  N <- nrow(mu)

  A <- matrix(NA, nrow = n, ncol = N)
  for (j in 1:N) {
    A[, j] <- purrr::map(
      1:n,
      ~ smoothing_kernel(as.numeric(x_obs[.x, ]), mu[j, ], sigma)
    ) %>%
      unlist()
  }
  A
}

#' Minimization Function for Rho Estimates
#'
#' @param x A numeric matrix of the unreduced dataset.
#' @param mu A numeric matrix of component centers.
#' @param w A numeric weight value.
#' @param lambda A given numeric vector.
#'
#' @return A numeric value.
#'
#' @noRd
f_lambda <- function(x, mu, w, lambda) {
  temp_denom <- Rfast::rowsums(t(t(cbind(rep(1, dim(mu)[1]), mu)) * lambda))
  temp_num <- mu * w

  f1 <- sum(w / temp_denom)
  f2 <- Rfast::colsums(temp_num / temp_denom)
  dist_euclidean(f1, 1) + dist_euclidean(f2, Rfast::colmeans(x))
}

#' Bound Vector Values Between 0 and 1
#'
#' @param theta A numeric vector.
#'
#' @return A numeric vector of bounded values.
#'
#' @noRd
bound_theta <- function(theta) {
  theta_bounded <- pmax(theta, 0)
  pmin(theta_bounded, 1)
}
