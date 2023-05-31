#' Principal Manifold Estimation
#'
#' This function still requires completed documentation
#'
#' @param x_obs Value
#' @param d Value
#' @param initialization Value
#' @param N0 Value
#' @param tuning_para_seq Value
#' @param alpha Value
#' @param max_comp Value
#' @param epsilon Value
#' @param max_iter Value
#' @param SSD_ratio_threshold Value
#' @param print_plots Value
#' @param verbose Value
#'
#' @return An object of type pme.
#' @export
#'
pme <- function(x_obs, d, initialization = NULL, N0 = 20 * D, tuning_para_seq = exp(-15:5), alpha = 0.05, max_comp = 100, epsilon = 0.05, max_iter = 100, SSD_ratio_threshold = 10, print_plots = FALSE, verbose = FALSE) {
  # Initial variable assignments --------------------------
  n <- dim(x_obs)[1]
  D <- dim(x_obs)[2]
  lambda <- 4 - d

  if (N0 == 0) {
    N0 <- 20 * D
  }
  if (N0 > max_comp) {
    max_comp <- 2 * N0
  }

  MSE_seq <- vector()
  SOL <- list()
  TNEW <- list()
  embeddings <- list()

  # Initialization ----------------------------------------

  if (is.null(initialization)) {
    initialization <- initialize_pme(x_obs, d, N0, alpha, max_comp)
  }
  initial_parameterization <- initialization$parameterization
  theta_hat <- initialization$theta_hat
  centers <- initialization$centers
  sigma <- initialization$sigma
  W <- diag(theta_hat)
  X <- centers
  I <- length(theta_hat)
  t_initial <- initial_parameterization$points

  # Fitting
  for (tuning_ind in 1:length(tuning_para_seq)) {
    w <- tuning_para_seq[tuning_ind]
    sol <- calc_coefficients(X, t_initial, W, w)

    fnew <- function(t) {
      as.vector(
        (t(sol[1:I, ]) %*% etaFunc(t, t_initial, lambda)) +
          (t(sol[(I + 1):(I + d + 1), ]) %*% matrix(c(1, t), ncol = 1))
      )
    }

    tnew <- calc_tnew(fnew, X, t_initial)
    SSD_new <- calc_SSD(fnew, X, tnew)

    count <- 1
    SSD_ratio <- 10 * epsilon

    while ((SSD_ratio > epsilon) && (SSD_ratio <= SSD_ratio_threshold) && (count <= max_iter)) {
      SSD_old <- SSD_new
      f0 <- fnew
      t_old <- tnew

      sol <- calc_coefficients(X, t_initial, W, w)

      fnew <- function(t) {
        as.vector(
          (t(sol[1:I, ]) %*% etaFunc(t, tnew, lambda)) +
            (t(sol[(I + 1):(I + d + 1), ]) %*% matrix(c(1, t), ncol = 1))
        )
      }

      tnew <- calc_tnew(fnew, X, t_old)
      SSD_new <- calc_SSD(fnew, X, tnew)

      SSD_ratio <- abs(SSD_new - SSD_old) / SSD_old
      count <- count + 1

      if (verbose == TRUE) {
        print_SSD(w, tuning_para_seq[tuning_ind], SSD_new, SSD_ratio, count)
      }
    }

    if (print_plots == TRUE) {
      plot_pme(fnew, x_obs, X, sol, tnew, I, d, lambda)
    }

    MSE_seq[tuning_ind] <- calc_msd(x_obs, initialization$km, fnew, tnew, D, d)
    SOL[[tuning_ind]] <- sol
    TNEW[[tuning_ind]] <- tnew
    embeddings[[tuning_ind]] <- fnew
    if (tuning_ind >= 4) {
      if (
        (MSE_seq[tuning_ind] > MSE_seq[tuning_ind - 1]) &&
        (MSE_seq[tuning_ind - 1] > MSE_seq[tuning_ind - 2]) &&
        (MSE_seq[tuning_ind - 2] > MSE_seq[tuning_ind - 3])
      ) {
        break
      }
    }
  }

  optimal_ind <- min(which(MSE_seq == min(MSE_seq)))
  sol_opt <- SOL[[optimal_ind]]
  t_opt <- TNEW[[optimal_ind]]

  if (verbose == TRUE) {
    paste0(
      "The optimal tuning parameter is ",
      as.character(tuning_para_seq[optimal_ind]),
      ", and the MSD of the optimal fit is ",
      as.character(MSE_seq[optimal_ind]),
      "."
    )
  }

  pme_out <- new_pme(
    embedding_map = embeddings[[optimal_ind]],
    knots = initialization$km,
    knot_weights = theta_hat,
    kernel_coefs = sol_opt[1:I, ],
    poly_coefs = sol_opt[(I + 1):(I + d + 1), ],
    tuning = tuning_para_seq[optimal_ind],
    MSD = MSE_seq,
    coefs = SOL,
    params = TNEW,
    tuning_vec = tuning_para_seq,
    embeddings = embeddings
  )
}

new_pme <- function(embedding_map,
                    knots,
                    knot_weights,
                    kernel_coefs,
                    poly_coefs,
                    tuning,
                    MSD,
                    coefs,
                    params,
                    tuning_vec,
                    embeddings) {
  pme_list <- list(
    embedding_map = embedding_map,
    knots = knots,
    knot_weights = knot_weights,
    kernel_coefs = kernel_coefs,
    poly_coefs = poly_coefs,
    tuning = tuning,
    MSD = MSD,
    coefs = coefs,
    params = params,
    tuning_vec = tuning_vec,
    embeddings = embeddings
  )
  vctrs::new_vctr(pme_list, class = "pme")
}

#' Is Object PME
#'
#' Check Whether Object Has Type PME
#'
#' @param x An object.
#'
#' @return Logical value.
#' @export
#'
#' @examples
#' num_value <- 5
#' is_pme(num_value)
is_pme <- function(x) {
  inherits(x, "pme")
}

initialize_pme <- function(x, d, N0, alpha, max_comp) {
  est <- hdmde(x, N0, alpha, max_comp)
  est_order <- order(est$mu[, 1])
  theta_hat <- est$theta_hat[est_order]
  centers <- est$mu[est_order, ]
  sigma <- est$sigma
  W <- diag(theta_hat)
  X <- est$mu[est_order, ]
  I <- length(theta_hat)

  dissimilarity <- as.matrix(stats::dist(X))
  init_parameterization <- vegan::isomap(dissimilarity, ndim = d, k = 10)

  list(
    parameterization = init_parameterization,
    theta_hat = theta_hat,
    centers = centers,
    sigma = sigma,
    km = est$km
  )
}

#' Compute Spline Coefficients
#'
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of low-dimensional paramterizations.
#' @param W Numeric matrix of cluster weights.
#' @param w smoothing parameter.
#'
#' @return
#' @export
#'
#' @examples
calc_coefficients <- function(X, t, W, w) {
  t_val <- cbind(rep(1, nrow(t)), t)
  E <- calcE(t, 4 - ncol(t))
  solve_weighted_spline(E, W, t_val, X, w, ncol(t), ncol(X))
}

#' Calculate a New Parameterization
#'
#' @param f Embedding map.
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of initial low-dimensional parameterizations.
#'
#' @return A numeric matrix of parameterizations.
#'
#' @noRd
calc_tnew <- function(f, X, t) {
  tnew <- purrr::map(1:nrow(X), ~projection_pme(X[.x, ], f, t[.x, ])) %>%
    unlist() %>%
    matrix(nrow = nrow(X), byrow = TRUE)
  tnew
}

#' Calculate Sum of Squared Distances
#'
#' @param f Embedding map.
#' @param X Numeric matrix of high-dimensional data.
#' @param t Numeric matrix of low-dimensional paramterizations.
#'
#' @return A numeric value.
#'
#' @noRd
calc_SSD <- function(f, X, t) {
  SSD_val <- purrr::map(1:nrow(X), ~dist_euclidean(X[.x, ], f(t[.x, ]))^2) %>%
    unlist() %>%
    sum()
  SSD_val
}

#' Print the Results of a PME Iteration
#'
#' @param w The smoothing parameter.
#' @param SSD_new The most recent SSD calculation.
#' @param SSD_ratio The ratio of the new and old SSD values.
#' @param count The iteration number.
#'
#' @noRd
print_SSD <- function(w, tuning_val, SSD_new, SSD_ratio, count) {
  print(
    paste0(
      "For tuning parameter ",
      as.character(round(tuning_val, 4)),
      ", iteration #",
      as.character(count),
      " gives SSD = ",
      as.character(round(SSD_new, 4)),
      " and SSD_ratio = ",
      as.character(round(SSD_ratio, 4)),
      "."
    )
  )
}

plot_pme <- function(f, x, centers, sol, t, I, d, lambda) {
  # pred_grid <- calc_tnew(centers, t, sol, I, d, lambda)
  pred_grid <- calc_tnew(f, centers, t)
  r_bounds <- Rfast::colMinsMaxs(pred_grid)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = nrow(centers)
    )
  }
  r_mat <- as.matrix(expand.grid(r_list))

  pred_grid <- r_mat
  f_pred <- purrr::map(
    1:nrow(pred_grid),
    ~f(unlist(as.vector(pred_grid[.x, ]))) %>%
      as.vector()
  ) %>%
    unlist() %>%
    matrix(nrow = nrow(pred_grid), byrow = TRUE)

  f_pred_full <- cbind(pred_grid, f_pred)

  if (dim(x)[2] == 2) {
    plt <- ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(
          x = x[, 1],
          y = x[, 2]
        ),
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = f_pred_full[, d + 1],
          y = f_pred_full[, d + 2]
        ),
        color = "red"
      )
    print(plt)
  } else if (dim(x)[2] >= 3) {
    plt <- plotly::plot_ly(
      x = f_pred_full[, d + 1],
      y = f_pred_full[, d + 2],
      z = f_pred_full[, d + 3],
      type = "scatter3d",
      mode = "markers",
      opacity = 0.5
    ) %>%
      plotly::add_markers(
        x = x[, 1],
        y = x[, 2],
        z = x[, 3],
        opacity = 0.15
      )
    print(plt)
  }
}

calc_msd <- function(x, km, f, t, D, d) {
  data_initial <- matrix(0, nrow = 1, ncol = D + d)
  for (i in 1:max(km$cluster)) {
    index_temp <- which(km$cluster == i)
    length_temp <- length(index_temp)
    temp_x <- matrix(x[index_temp, ], nrow = length_temp)
    t_temp <- matrix(rep(t[i, 1], length_temp))
    for (j in 1:d) {
      t_temp <- cbind(t_temp, rep(t[i, j], length_temp))
    }
    t_temp <- matrix(t_temp[, -1], nrow = length_temp)
    data_initial <- rbind(data_initial, cbind(temp_x, t_temp))
  }
  data_initial <- data_initial[-1, ]
  proj_para <- purrr::map(
    1:nrow(data_initial),
    ~projection_pme(data_initial[.x, 1:D], f, data_initial[.x, (D + 1):(D + d)]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)
  proj_points <- purrr::map(
    1:nrow(proj_para),
    ~f(proj_para[.x, ]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)

  mse <- purrr::map(
    1:nrow(data_initial),
    ~dist_euclidean(data_initial[.x, 1:D], proj_points[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  mse
}
