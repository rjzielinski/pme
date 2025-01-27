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

#' High-Dimensional Mixture Density Estimation - Modification
#'
#' Data reduction approach to estimate a high-dimensional dataset using a
#' mixture of distributions. The centers of these distributions are identified
#' using k-means clustering. This function is modified to use the BIC of the
#' k-means clustering output to select the number of mixture components.
#'
#' @param x_obs A numeric matrix containing the data to be reduced.
#' @param N0 An integer specifying the lower bound for the number of density components.
#' @param alpha A numeric value between 0 and 1 specifying the confidence level.
#' @param max_comp An integer specifying the upper bound for the number of density components.
#'
#' @return An object of type hdmde
#' @export
hdmde_mod <- function(x_obs, N0, alpha, max_comp) {
  # Initialization ----------------------------
  n <- nrow(x_obs)
  D <- ncol(x_obs)
  N <- N0

  component_estimates <- compute_estimates(x_obs, N)
  km_old <- component_estimates$km

  # see http://sherrytowers.com/2013/10/24/k-means-clustering/ for aic
  # computations

  aic_old <- km_old$tot.withinss + (2 * N * D)
  test_rejection <- FALSE

  while ((test_rejection == FALSE) & (N < min(n - 1, max_comp))) {
    N <- N + 1
    components_new <- compute_estimates(x_obs, N)
    km_new <- components_new$km
    aic_new <- km_new$tot.withinss + (2 * N * D)

    if (abs(aic_new - aic_old) / aic_old < alpha) {
      test_rejection <- TRUE
    }
    km_old <- km_new
    aic_old <- aic_new
    component_estimates <- components_new
  }

  new_hdmde(components_new$mu, components_new$sigma, components_new$theta_hat, components_new$km)
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

  log_p_old <- purrr::map(
    1:n,
    ~ f_test(
      x_obs[.x, ],
      component_estimates$mu,
      component_estimates$sigma,
      component_estimates$theta_hat
    )
  ) %>%
    unlist()

  test_rejection <- FALSE

  while ((test_rejection == FALSE) & (N < min(n - 1, max_comp))) {
    N <- N + 1
    components_new <- compute_estimates(x_obs, N)
    log_p_new <- purrr::map(
      1:n,
      ~ f_test(
        x_obs[.x, ],
        components_new$mu,
        components_new$sigma,
        components_new$theta_hat
      )
    ) %>%
      unlist()

    log_difference <- purrr::map(
      1:n,
      ~ logspace_diff(log_p_new[.x], log_p_old[.x])
    ) %>%
      unlist()

    diff_ind <- log_p_new > log_p_old

    log_diff_sq <- purrr::map(
      1:n,
      ~ logspace_diff(logspace_sum(2 * log_p_new[.x], 2 * log_p_old[.x]), log(2) + log_p_new[.x] + log_p_old[.x])
    ) %>%
      unlist()

    log_mean_diff <- 0
    p_new_improve <- which(diff_ind)
    p_old_improve <- which(!diff_ind)
    log_mean_diff <- logspace_diff(
      logspace_sum_vec(log_difference[p_new_improve]),
      logspace_sum_vec(log_difference[p_old_improve])
    ) - log(n)

    # we have log difference values, but they are 0 after exponentiating
    # difference <- p_new - p_old

    # log_mean_diff <- logspace_sum_vec(log_difference) - log(length(log_difference))
    log_residuals <- purrr::map(1:n, ~ logspace_diff(log_difference[.x], log_mean_diff)) %>%
      unlist()
    # log_sigma_hat_sq <- 2 * logspace_sum_vec(log_difference) - 2 * logspace_sum_vec(log_mean_diff)
    log_sigma_hat_sq <- purrr::map(
      1:n,
      ~ logspace_diff(log_diff_sq[.x], 2 * log_mean_diff)
    ) %>%
      unlist() %>%
      logspace_sum_vec()
    log_sigma_hat_sq <- log_sigma_hat_sq - log(n)

    # sigma_hat_sq <- mean((difference - mean(difference))^2)
    # Z_I_N <- sqrt(n) * mean(difference) / sqrt(sigma_hat_sq)
    Z_I_N <- exp(log_mean_diff + (0.5 * log(n)) - (0.5 * log_sigma_hat_sq))
    if (!is.na(Z_I_N) & (abs(Z_I_N) <= zalpha)) {
      test_rejection <- TRUE
    }
    log_p_old <- log_p_new
  }

  new_hdmde(components_new$mu, components_new$sigma, components_new$theta_hat, components_new$km)
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
  km <- stats::kmeans(
    x,
    n,
    iter.max = 10000,
    nstart = 1000,
    algorithm = "MacQueen"
  )
  mu <- km$centers
  sigma_est <- estimate_sigma(x, km)
  # theta_hat <- calc_weights(x, mu, sigma_est)
  # instead of using calc_weights function, approximate weights using
  # kmeans cluster assignments for stability in high dimensions
  theta_hat <- purrr::map(
    1:n,
    ~ sum(km$cluster == .x) / nrow(x)
  ) %>%
    unlist()

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
  clusters <- unique(km$cluster)
  sigma_vec <- rep(NA, length(clusters))
  for (j in seq_along(clusters)) {
    index_temp <- which(km$cluster == clusters[j])
    x_temp <- matrix(x[index_temp, ], nrow = length(index_temp))
    s <- purrr::map(
      seq_len(nrow(x_temp)),
      ~ dist_euclidean(x_temp[.x, ], km$centers[clusters[j], ])^2
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
  A <- calc_A(x_obs, mu, sigma)
  theta_old <- rep(1 / N, N)
  abs_diff <- 10 * epsilon
  count <- 0
  lambda_hat_old <- c(n, rep(-1, D))


  while ((abs_diff > epsilon) & (count <= max_iter)) {
    W <- t(t(A) + log(theta_old))
    W_adj <- W - apply(W, 1, logspace_sum_vec)
    w <- apply(W_adj, 2, logspace_sum_vec)
    opt_out <- stats::optim(
      par = lambda_hat_old,
      fn = f_lambda,
      x = x_obs,
      mu = mu,
      w = exp(w),
      method = "CG",
      control = list(maxit = 1000)
    )
    lambda_hat2 <- opt_out$par
    lambda_hat <- stats::nlm(
      f = f_lambda,
      p = lambda_hat_old,
      x = x_obs,
      mu = mu,
      w = exp(w),
      iterlim = 10000
    )$estimate
    theta_new <- exp(w) / Rfast::rowsums(t(t(cbind(rep(1, N), mu)) * lambda_hat))
    theta_new2 <- theta_new / sum(theta_new)

    abs_diff <- dist_euclidean(theta_new, theta_old)
    # should I change how theta is bounded?
    # does theta need to be normalized?
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
    1:nrow(mu),
    ~ log_smoothing_kernel(x, mu[.x, ], sigma)
  ) %>%
    unlist()
  logspace_sum_vec(log(theta_hat) + comp_vec)
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
