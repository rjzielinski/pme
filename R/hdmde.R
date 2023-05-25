#' Create New HDMDE Object
#'
#' @param mu A numeric matrix containing the density component centers.
#' @param sigma A numeric value describing the cluster variance.
#' @param theta A numeric vector of weight values.
#'
#' @return An object of type hdmde.
#'
#' @noRd
new_hdmde <- function(mu, sigma, theta) {
  hdmde_list <- list(
    mu = mu,
    sigma = sigma,
    theta_hat = theta_hat
  )
  vctrs::new_vctr(hdmde_list, class = "hdmde")
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
#'
#' @examples
hdmde <- function(x_obs, N0, alpha, max_comp) {
  # Initialization ----------------------------
  zalpha <- qnorm(1 - (alpha / 2))
  n <- nrow(x_obs)
  D <- ncol(x_obs)
  N <- N0

  component_estimates <- compute_estimates(x_obs, N)
  p_old <- map(
    1:n,
    ~ f_test(
      x_obs[.x, ],
      component_estimates$mu,
      component_estimates$sigma,
      component_estimates$theta_hat
    )
  )

  test_rejection <- FALSE

  while ((test_rejection == FALSE) & (N <= min(n, max_comp))) {
    N <- N + 1
    components_new <- compute_estimates(x_obs, N)
    p_new <- map(
      1:n,
      ~ f_test(
        x_obs[.x, ],
        components_new$mu,
        components_new$sigma,
        components_new$theta_hat
      )
    )
    difference <- p_new - p_old
    sigma_hat_sq <- mean((difference - mean(difference))^2)
    Z_I_N <- sqrt(n) * mean(difference) / sqrt(sigma_hat_sq)
    if (!is.na(Z_I_N) & (abs(Z_I_N) <= zalpha)) {
      test_rejection <- TRUE
    }
    p_old <- p_new
  }


}


# HDMDE HELPER FUNCTIONS -------------------------------------------------------


#' Compute Mixture Component Estimates
#'
#' For a given value of N, estimate the mixture component centers, variance, and weights.
#'
#' @param x_obs A numeric matrix.
#' @param N The number of mixture components to be estimated.
#'
#' @return A list containing the mixture centers, variance, and weights.
#'
#' @noRd
compute_estimates <- function(x_obs, N) {
  # increasing iter.max helps to avoid poor fit in high dimensional situations.
  km <- kmeans(x_obs, N, iter.max = 10000, nstart = 100)
  mu <- km$centers
  sigma_est <- estimate_sigma(km)
  theta_hat <- weight_seq(x_obs, mu, sigma_est)

  list(
    mu = mu,
    sigma = sigma_est,
    theta_hat = theta_hat
  )
}

#' Estimate Mixture Component Variance
#'
#' @param km A k-means object including the estimated mixture components.
#'
#' @return A numeric estimate of the variance of the mixture components.
#'
#' @noRd
estimate_sigma <- function(km) {
  N <- nrow(km$centers)
  sigma_vec <- rep(NA, N)
  for (j in 1:N) {
    index_temp <- which(km$cluster == j)
    x_temp <- x_obs[index_temp, ]
    s <- map(
      1:nrow(x_temp),
      ~ dist_euclidean(x_temp[.x, ], mu[j, ])^2
    ) %>%
      unlist()
    sigma_vec[j] <- mean(s)
  }
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
  theta_old <- rep(1 / N, N)
  abs_diff <- 10 * epsilon
  count <- 0
  lambda_hat_old <- c(n, rep(-1, D))

  A <- calc_A(x_obs, mu, sigma)

  while ((abs_diff > epsilon) & (count <= max_iter)) {
    W <- t(t(A) * theta_old)
    w <- Rfast::colsums(W / Rfast::rowsums(W))
    lambda_hat <- nlm(f_lambda, lambda_hat_old, iterlim = 1000)$estimate
    theta_new <- w / Rfast::rowsums(t(t(cbind(rep(1, N), mu)) * lambda_hat))

    abs_diff <- dist_euclidean(theta_new, theta_old)
    if (is.na(abs_diff)) {
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
  comp_vec <- map(
    1:N,
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
calc_A <- function(x_obs, mu, sigma) {
  n <- nrow(x_obs)
  N <- nrow(mu)

  A <- matrix(NA, nrow = n, ncol = N)
  for (j in 1:N) {
    A[, j] <- map(
      1:n,
      ~ ker(x_obs[.x, ], mu[j, ], sigma)
    ) %>%
      unlist()
  }
  A
}

#' Minimization Function for Rho Estimates
#'
#' @param mu A numeric matrix of component centers.
#' @param w A numeric weight value.
#' @param lambda A given numeric vector.
#'
#' @return A numeric value.
#'
#' @noRd
f_lambda <- function(mu, w, lambda) {
  temp_denom <- Rfast::rowsums(t(t(cbind(rep(1, dim(mu)[1])) * lambda)))
  temp_num <- mu * w

  f1 <- sum(w / temp_denom)
  f2 <- Rfast::colsums(temp_num / temp_denom)
  dist_euclidean(f1, 1) + dist_euclidean(f2, Rfast::colmeans(x_obs))
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
