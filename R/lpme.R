#' Create New LPME Object
#'
#' @param embedding_map The optimal embedding function.
#' @param msd Vector of Mean Squared Distance Values.
#' @param coefficients A list of matrices of spline coefficients.
#' @param times Numeric vector of the time values included in the dataset.
#' @param initial_parameterization A matrix of the initial parameters found for the data.
#' @param optimal_parameterization A matrix representing the best set of parameters.
#' @param d Integer value representing the intrinsic dimension of the manifold.
#' @param D Integer value representing the dimension of the input data.
#' @param n_knots The number of components identified.
#' @param lambda A vector of smoothing parameters for the pme functions at each time point.
#' @param gamma A vector of smoothing parameter values for the higher level spline model.
#' @param coefficient_list A list containing all matrices of coefficients.
#' @param coefficient_functions A list of functions that produce the PME coefficients at a given time point.
#' @param parameterization_list A list containing all matrices of parameters.
#' @param smoothing_method Choose between smoothing by "spline" or "gp".
#' @param initialization_algorithm The algorithm that was used to identify the initial parameter values.
#'
#' @return An object with class "lpme".
#' @export
new_lpme <- function(
  embedding_map,
  msd,
  coefficients,
  times,
  initial_parameterization,
  optimal_parameterization,
  d,
  D,
  n_knots,
  lambda,
  gamma,
  coefficient_list,
  coefficient_functions,
  parameterization_list,
  smoothing_method,
  initialization_algorithm
) {
  lpme_list <- list(
    embedding_map = embedding_map,
    d = d,
    D = D,
    lambda = lambda,
    gamma = gamma,
    times = times,
    n_knots = n_knots,
    initial_parameterization = initial_parameterization,
    optimal_parameterization = optimal_parameterization,
    sol_coef = coefficients,
    msd = msd,
    sol_coef_list = coefficient_list,
    sol_coef_functions = coefficient_functions,
    parameterization_list = parameterization_list,
    smoothing_method = smoothing_method,
    initialization = initialization_algorithm
  )
  vctrs::new_vctr(lpme_list, class = "lpme")
}

#' Check Whether Object is of Type LPME
#'
#' @param x An object to test.
#'
#' @return A logical value.
#' @export
is_lpme <- function(x) {
  inherits(x, "lpme")
}

#' Fit an LPME Object
#'
#' @param data A numeric matrix of data.
#' @param d The intrinsic dimension of the data.
#' @param smoothing_method The approach taken to smoothing over PME coefficients.
#' @param gamma A vector of numeric smoothing parameter values.
#' @param lambda A vector of numeric smoothing parameter values for the PME algorithm.
#' @param initialization_algorithm A character value describing the algorithm used to obtain initial parameter values. Accepted values are "isomap", "diffusion_maps", and "laplacian_eigenmaps".
#' @param init_type A character value describing whether initialization should use cluster centers or subsample from cluster values. Accepted values are "centers" or "subsampling".
#' @param alpha A value.
#' @param min_clusters The minimum number of clusters in the data.
#' @param max_clusters The maximum number of clusters identified in the data.
#' @param epsilon A value.
#' @param max_iter The maximum number of iterations.
#' @param verbose A logical value indicating whether messages should be printed.
#' @param print_plots A logical value indicating whether plots should be printed.
#' @param increase_threshold A value.
#' @param init Indicates which time points are used to initialize the function.
#'
#' @return An object of type "lpme".
#' @export
lpme <- function(
  data,
  d,
  smoothing_method = "spline",
  gamma = NULL,
  lambda = NULL,
  initialization_algorithm = "isomap",
  init_type = "centers",
  alpha = 0.05,
  min_clusters = 0,
  max_clusters = 500,
  epsilon = 0.05,
  max_iter = 100,
  verbose = TRUE,
  print_plots = TRUE,
  increase_threshold = 1.05,
  init = "full"
) {
  # Declare initial variable values ---------------------------------------
  time_points <- unique(data[, 1])
  min_observations <- nrow(data)
  for (time_val in time_points) {
    temp_data <- data[data[, 1] == time_val, ]
    min_observations <- min(min_observations, nrow(temp_data))
  }

  max_clusters <- min(min_observations - 1, max_clusters)

  if (is.null(gamma)) {
    if (smoothing_method == "spline") {
      gamma <- c(0, exp(-15:10))
    } else if (smoothing_method == "gp") {
      gamma <- 0:10 + 0.5
    }
  }

  if (is.null(lambda)) {
    lambda <- exp(-20:10)
  }

  if (min_clusters == 0) {
    min_clusters <- 10 * d
  }

  initialization <- initialize_lpme(
    data,
    init,
    time_points,
    d,
    alpha,
    max_clusters,
    min_clusters,
    initialization = initialization_algorithm
  )

  init_pme_list <- fit_init_pmes(
    data,
    time_points,
    init,
    initialization,
    d,
    lambda
  )
  splines <- merge_spline_coefs(init_pme_list, d, time_points)

  spline_coefficients <- splines$coef_full
  x_merged <- splines$x_test
  params <- splines$params
  times <- splines$times
  n_knots <- splines$n_knots
  lambda <- splines$lambda

  D_coef <- dim(spline_coefficients)[2]
  D_out <- dim(x_merged)[2]
  d <- dim(params)[2]
  n <- dim(x_merged)[1]
  I_new <- n
  t_initial <- params %>%
    as.matrix()

  MSE_seq_new <- vector()
  SOL_coef <- list()
  TNEW_new <- list()
  coefs <- list()
  x_funs <- list()
  functions <- list()
  func_coef <- list()

  inv_errors <- 1 / init_pme_list$errors
  weights <- inv_errors / sum(inv_errors)

  for (tuning_ind in 1:length(gamma)) {
    if (smoothing_method == "spline") {
      f_coef_list <- compute_f_coef(
        gamma[tuning_ind],
        diag(weights),
        t_initial,
        times,
        spline_coefficients,
        1,
        D_coef,
        4 - d,
        3
      )

      f_new <- function(t) {
        coefs <- f_coef_list$f(t[1])
        coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
        return_vec <- t(coef_mat[1:n_knots, ]) %*%
          etaFunc(t[-1], t_initial, 4 - d) +
          t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*%
            matrix(c(1, t[-1]), ncol = 1)
        c(t[1], return_vec)
      }
    } else if (smoothing_method == "gp") {
      # invisible(
      #   capture.output({
      #     gp <- GPFDA::gpr(
      #       response = spline_coefficients,
      #       input = times,
      #       Cov = "matern",
      #       meanModel = "t",
      #       nu = tuning_para_seq[tuning_ind]
      #     )
      #   })
      # )
      #
      # f_new <- function(t) {
      #   coefs <- GPFDA::gprPredict(
      #     train = gp,
      #     inputNew = t[1],
      #     noiseFreePred = TRUE
      #   )$pred.mean %>%
      #     as.vector()
      #   coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
      #   return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, 4 - d) +
      #     t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
      #   c(t[1], return_vec)
      # }
    }

    updated_param <- update_parameterization(
      time_points,
      t_initial,
      x_merged,
      f_new,
      n_knots,
      d,
      d,
      4 - d
    )

    if (print_plots == TRUE) {
      plot_lpme(data, f_new, t_initial, d, D_out - 1, time_points)
    }

    # Cross-validation section
    nearest_x <- calc_nearest_x(data, x_merged)
    init_param <- calc_init_param(
      data,
      updated_param$parameterization,
      nearest_x
    )

    cv_mse <- calc_mse_cv(
      leave_one_out = TRUE,
      f = f_new,
      df = data,
      init_param = init_param,
      time_points = time_points,
      r = updated_param$parameterization[, -1],
      r_initial = t_initial,
      n_knots = n_knots,
      d = d,
      d_new2 = 1,
      D_out = D_out - 1,
      D_coef = D_coef,
      lambda = lambda,
      gamma = 4 - d,
      gamma2 = 3,
      r_full2 = times,
      w = gamma[tuning_ind],
      smoothing_method = smoothing_method
    )

    data_n <- sapply(time_points, function(x) nrow(data[data[, 1] == x, ]))
    MSE_new <- stats::weighted.mean(cv_mse, data_n)
    MSE_seq_new[tuning_ind] <- MSE_new

    if (verbose == TRUE) {
      print_mse(gamma[tuning_ind], MSE_new)
    }

    if (smoothing_method == "spline") {
      SOL_coef[[tuning_ind]] <- f_coef_list$sol
    } else {
      SOL_coef[[tuning_ind]] <- NA
    }

    TNEW_new[[tuning_ind]] <- updated_param$parameterization
    coefs[[tuning_ind]] <- params
    x_funs[[tuning_ind]] <- updated_param$embedding
    functions[[tuning_ind]] <- f_new
    func_coef[[tuning_ind]] <- ifelse(
      smoothing_method == "spline",
      f_coef_list$f,
      NA
    )

    if (tuning_ind >= 8) {
      if (!is.unsorted(MSE_seq_new[(tuning_ind - 5):tuning_ind])) {
        break
      }
    }
  }

  optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
  sol_opt <- SOL_coef[[optimal_ind]]
  t_new_opt <- TNEW_new[[optimal_ind]]
  coefs_opt <- coefs[[optimal_ind]]
  f.optimal <- functions[[optimal_ind]]

  if (verbose == TRUE) {
    plot_MSE(MSE_seq_new, gamma, optimal_ind)
  }

  lpme_out <- new_lpme(
    embedding_map = f.optimal,
    msd = MSE_seq_new,
    coefficients = sol_opt,
    times = as.vector(time_points),
    initial_parameterization = t_initial,
    optimal_parameterization = t_new_opt,
    d = d,
    D = D_out,
    n_knots = n_knots,
    lambda = lambda,
    gamma = gamma,
    coefficient_list = SOL_coef,
    coefficient_functions = func_coef,
    parameterization_list = TNEW_new,
    smoothing_method = smoothing_method,
    initialization_algorithm = initialization_algorithm
  )

  lpme_out
}

#' Embed Low-Dimensional Parameterization in High-Dimensional Space
#'
#' @param object An object of class `lpme`.
#' @param x A value
#'
#' @return A value
#' @export
#'
embed <- function(object, x) {
  n_times <- length(object$times)
  nu <- 4 - object$d
  f_coef <- function(t) {
    vec <- as.vector(
      t(object$sol_coef[1:n_times, ]) %*%
        etaFunc(t, as.matrix(object$times, ncol = 1), 3) +
        t(object$sol_coef[(n_times + 1):(n_times + 2), ]) %*%
          matrix(c(1, t), ncol = 1)
    )
    return(vec)
  }

  coefs <- f_coef(x[1])
  coef_mat <- matrix(coefs, object$n_knots + object$d + 1, byrow = TRUE)
  vec <- as.vector(
    t(coef_mat[1:object$n_knots, ]) %*%
      etaFunc(x[-1], object$r_init, nu) +
      t(coef_mat[(object$n_knots + 1):(object$n_knots + object$d + 1), ]) %*%
        matrix(c(1, x[-1]), ncol = 1)
  )
  return(c(x[1], vec))
}

parameterize <- function(object, x) {
  embeddings <- purrr::map(
    1:nrow(object$r_fit),
    ~ embed(object, object$r_fit[.x, ])
  ) %>%
    purrr::reduce(rbind)
  nearest <- calc_nearest_x(
    matrix(x, nrow = 1),
    embeddings[embeddings[, 1] == x[1], ]
  )
  estimate <- projection_lpme(
    x,
    function(t) embed(object, t),
    object$r_fit[object$r_fit[, 1] == x[1], ][nearest, ]
  )

  return(estimate)
}

initialize_lpme <- function(
  df,
  init,
  time_points,
  d,
  alpha,
  max_comp,
  min_comp,
  initialization,
  initialization_type
) {
  if (init %in% c("first", "full")) {
    if (init == "first") {
      init_df <- df[df[, 1] == time_points[1], -1]
      init_pme <- pme(
        init_df,
        d,
        initialization_algorithm = initialization,
        initialization_type = initialization_type
      )
      init_pme_center_order <- order(init_pme$knots$centers[, 1])
      init_pme_centers <- init_pme$knots$centers[init_pme_center_order, ]

      init_timevals <- list()
      init_theta_hat <- list()
      init_centers <- list()
      init_sigma <- list()
      init_clusters <- list()
      init_parameterization <- list()
      init_D <- ncol(df[, -1])
      init_n <- nrow(df)
      lambda <- 4 - d
      init_N0 <- min_comp

      for (idx in 1:length(time_points)) {
        init_df_temp <- df[df[, 1] == time_points[idx], -1]
        init_est_temp <- hdmde(init_df_temp, init_N0, alpha, max_comp)
        est_temp_order <- order(init_est_temp$mu[, 1])
        init_timevals[[idx]] <- rep(time_points[idx], dim(init_est_temp$mu)[1])
        init_theta_hat[[idx]] <- init_est_temp$theta_hat[est_temp_order]
        init_centers[[idx]] <- init_est_temp$mu[est_temp_order, ]
        init_sigma[[idx]] <- init_est_temp$sigma
        init_clusters[[idx]] <- init_est_temp$km

        nearest_clusters <- purrr::map(
          1:nrow(init_centers[[idx]]),
          ~ which.min(as.vector(apply(
            init_pme_centers,
            1,
            dist_euclidean,
            y = init_centers[[idx]][.x, ]
          )))
        ) %>%
          purrr::reduce(c)

        opt_run <- which.min(init_pme$MSD)
        params_init <- init_pme$parameterization[[opt_run]][
          nearest_clusters,
        ] %>%
          matrix(nrow = nrow(init_centers[[idx]]), byrow = TRUE)

        init_parameterization[[idx]] <- purrr::map(
          1:nrow(init_centers[[idx]]),
          ~ projection_pme(
            init_centers[[idx]][.x, ],
            init_pme$embedding_map,
            params_init[.x, ]
          )
        ) %>%
          purrr::reduce(rbind)
      }

      init_list <- list(
        timevals = init_timevals,
        theta_hat = init_theta_hat,
        centers = init_centers,
        sigma = init_sigma,
        clusters = init_clusters,
        isomap = init_parameterization
      )
    } else if (init == "full") {
      init_timevals <- list()
      init_theta_hat <- list()
      init_centers <- list()
      init_sigma <- list()
      init_clusters <- list()

      init_dimension_size <- dim(df[, -1])
      init_D <- init_dimension_size[2]
      init_n <- init_dimension_size[1]
      lambda <- 4 - d
      if (is.null(min_comp)) {
        init_N0 <- 10 * init_D
      } else {
        init_N0 <- min_comp
      }
      for (idx in 1:length(time_points)) {
        init_df_temp <- df[df[, 1] == time_points[idx], -1]
        init_est_temp <- hdmde(init_df_temp, init_N0, alpha, max_comp)
        est_temp_order <- order(init_est_temp$mu[, 1])
        init_timevals[[idx]] <- rep(time_points[idx], dim(init_est_temp$mu)[1])
        init_theta_hat[[idx]] <- init_est_temp$theta_hat[est_temp_order]
        init_centers[[idx]] <- init_est_temp$mu[est_temp_order, ]
        init_sigma[[idx]] <- init_est_temp$sigma
        init_clusters[[idx]] <- init_est_temp$km
      }

      init_timevals <- purrr::reduce(init_timevals, c)
      init_theta_hat <- purrr::reduce(init_theta_hat, c)
      init_centers <- purrr::reduce(init_centers, rbind)
      init_sigma <- purrr::reduce(init_sigma, c)
      init_W <- diag(init_theta_hat)
      init_X <- init_centers
      init_I <- length(init_theta_hat)

      init_dissimilarity_matrix <- as.matrix(stats::dist(init_X))
      init_isomap <- vegan::isomap(
        init_dissimilarity_matrix,
        ndim = d,
        k = floor(sqrt(nrow(init_dissimilarity_matrix)))
      )
      init_list <- list(
        timevals = init_timevals,
        theta_hat = init_theta_hat,
        centers = init_centers,
        sigma = init_sigma,
        clusters = init_clusters,
        isomap = init_isomap
      )
    }

    return(init_list)
  } else {
    print("Please choose `init` option 'first' or 'full'.")
  }
}

fit_init_pmes <- function(df, time_points, init_option, init, d, lambda) {
  pme_results <- list()
  kernel_coefs <- list()
  polynomial_coefs <- list()
  funcs <- list()
  clusters <- list()
  centers <- list()
  embeddings <- list()
  params <- list()
  times <- list()
  num_clusters <- rep(0, length(time_points))
  errors <- vector()

  for (idx in seq_along(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    if (init_option == "full") {
      pme_results[[idx]] <- pme(
        data = df_temp[, -1],
        d = d,
        lambda = lambda,
        initialization = list(
          parameterization = matrix(
            init$isomap$points[init$timevals == time_points[idx], ],
            nrow = length(init$theta_hat[init$timevals == time_points[idx]])
          ),
          theta_hat = init$theta_hat[init$timevals == time_points[idx]],
          centers = init$centers[init$timevals == time_points[idx], ],
          sigma = init$sigma[idx],
          km = init$clusters[[idx]]
        ),
        verbose = FALSE,
        print_plots = FALSE
      )
    } else if (init_option == "first") {
      pme_results[[idx]] <- pme(
        data = df_temp[, -1],
        d = d,
        lambda = lambda,
        initialization = list(
          parameterization = init$isomap[[idx]],
          theta_hat = init$theta_hat[[idx]],
          centers = init$centers[[idx]],
          sigma = init$sigma[[idx]],
          km = init$clusters[[idx]]
        ),
        verbose = FALSE,
        print_plots = FALSE
      )
    } else {
      pme_results[[idx]] <- pme(
        data = df_temp[, -1],
        d = d,
        lambda = lambda,
        verbose = FALSE
      )
    }

    opt_run <- which.min(pme_results[[idx]]$MSD)
    funcs[[idx]] <- pme_results[[idx]]$embedding_map
    centers[[idx]] <- pme_results[[idx]]$knots$centers
    num_clusters[idx] <- dim(pme_results[[idx]]$knots$centers)[1]
    if (idx == 1) {
      clusters[[idx]] <- pme_results[[idx]]$knots$cluster
    } else {
      clusters[[idx]] <- pme_results[[idx]]$knots$cluster + sum(num_clusters)
    }
    kernel_coefs[[idx]] <- pme_results[[idx]]$coefs[[opt_run]][
      1:num_clusters[idx],
    ]
    polynomial_coefs[[idx]] <- pme_results[[idx]]$coefs[[opt_run]][
      (num_clusters[idx] + 1):(num_clusters[idx] + d + 2),
    ] %>%
      t() %>%
      matrix(nrow = 1)
    embeddings[[idx]] <- apply(
      pme_results[[idx]]$parameterization[[opt_run]],
      1,
      funcs[[idx]]
    ) %>%
      t()
    params[[idx]] <- pme_results[[idx]]$parameterization[[opt_run]]
    times[[idx]] <- time_points[idx]
    errors[idx] <- pme_results[[idx]]$MSD[opt_run]
  }

  init_pme_list <- list(
    pme_results = pme_results,
    kernel_coefs = kernel_coefs,
    polynomial_coefs = polynomial_coefs,
    funcs = funcs,
    clusters = clusters,
    centers = centers,
    embeddings = embeddings,
    params = params,
    times = times,
    num_clusters = num_clusters,
    errors = errors
  )
  init_pme_list
}

merge_spline_coefs <- function(pme_list, d, time_points) {
  lambda <- vector()
  x_vals <- list()
  spline_coefs <- list()

  params <- purrr::reduce(pme_list$params, rbind)
  length_param <- ceiling(max(pme_list$num_clusters)^(1 / d))
  params <- gen_parameterization(params, ceiling(max(pme_list$num_clusters)), d)
  n_knots <- nrow(params)
  times <- purrr::reduce(pme_list$times, rbind)

  for (time_idx in 1:length(time_points)) {
    lambda[time_idx] <- pme_list$pme_results[[time_idx]]$tuning
    f <- pme_list$funcs[[time_idx]]
    output <- purrr::map(
      1:n_knots,
      ~ f(params[.x, ])
    ) %>%
      purrr::reduce(rbind)
    x_vals[[time_idx]] <- cbind(time_points[time_idx], output)

    R <- cbind(rep(1, n_knots), params)
    E <- calcE(params, 4 - d)

    spline_coefs[[time_idx]] <- solve_spline(
      E,
      R,
      output,
      lambda[time_idx],
      ncol(params),
      ncol(output)
    ) %>%
      t() %>%
      matrix(nrow = 1)
  }
  coef_full <- purrr::reduce(spline_coefs, rbind)
  x_test <- purrr::reduce(x_vals, rbind)

  out_list <- list(
    params = params,
    times = times,
    coef_full = coef_full,
    x_test = x_test,
    n_knots = n_knots,
    lambda = lambda
  )
  out_list
}

compute_f_coef <- function(
  w_new,
  W,
  t_new,
  r_full2,
  coef_full,
  d_new2,
  D_coef,
  gamma,
  gamma2
) {
  T_new <- cbind(rep(1, nrow(t_new)), t_new)
  T_new2 <- cbind(rep(1, nrow(r_full2)), r_full2)

  E_new <- calcE(t_new, gamma)
  E_new2 <- calcE(r_full2, gamma2)

  sol_coef <- solve_weighted_spline(
    E_new2,
    W,
    T_new2,
    coef_full,
    w_new,
    d_new2,
    D_coef
  )
  f_coef <- function(t) {
    return_vec <- as.vector(
      t(sol_coef[1:nrow(coef_full), ]) %*%
        etaFunc(t, r_full2, gamma2) +
        t(sol_coef[(nrow(coef_full) + 1):(nrow(coef_full) + d_new2 + 1), ]) %*%
          matrix(c(1, t), ncol = 1)
    )
    return_vec
  }
  out_list <- list(
    f = f_coef,
    sol = sol_coef
  )
}

plot_lpme <- function(x, f, r, d, D, time_points) {
  time_vals <- seq(min(time_points), max(time_points), 0.05)
  new_pred <- matrix(ncol = ncol(r) + 1)
  for (time in time_vals) {
    new_pred <- rbind(new_pred, cbind(time, r))
  }
  pred_grid <- new_pred[-1, ]

  f_pred <- purrr::map(
    1:nrow(pred_grid),
    ~ f(pred_grid[.x, ])
  ) %>%
    purrr::reduce(rbind)
  f_pred_full <- cbind(pred_grid, f_pred)

  if (D == 2) {
    plt <- plotly::plot_ly(
      x = f_pred_full[, d + 3],
      y = f_pred_full[, d + 4],
      z = f_pred_full[, 1],
      type = "scatter3d",
      mode = "markers",
      marker = list(size = 3)
    ) %>%
      plotly::add_markers(
        x = x[, 2],
        y = x[, 3],
        z = x[, 1],
        opacity = 0.15
      )
    print(plt)
  } else {
    plt <- plotly::plot_ly(
      x = f_pred_full[, d + 3],
      y = f_pred_full[, d + 4],
      z = f_pred_full[, d + 5],
      frame = f_pred_full[, d + 2],
      type = "scatter3d",
      mode = "markers",
      opacity = 1,
      marker = list(size = 3)
    ) %>%
      plotly::add_markers(
        x = x[, 2],
        y = x[, 3],
        z = x[, 4],
        frame = x[, 1],
        opacity = 0.2
      ) %>%
      plotly::layout(
        scene = list(
          xaxis = list(range = list(min(x[, 2]), max(x[, 2]))),
          yaxis = list(range = list(min(x[, 3]), max(x[, 3]))),
          zaxis = list(range = list(min(x[, 4]), max(x[, 4])))
        )
      )
    print(plt)
  }
}

plot_MSE <- function(MSE_seq_new, gamma, optimal_ind) {
  plt <- ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(
        x = log(gamma[1:length(MSE_seq_new)]),
        y = MSE_seq_new
      )
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = log(gamma[1:length(MSE_seq_new)]),
        y = MSE_seq_new
      ),
      color = "orange"
    ) +
    ggplot2::geom_vline(
      xintercept = log(gamma[optimal_ind]),
      color = "darkgreen"
    ) +
    ggplot2::xlab("Log Gamma") +
    ggplot2::ylab("MSD")
  print(
    paste0(
      "The optimal tuning parameter is ",
      as.character(gamma[optimal_ind]),
      ", and the MSD of the optimal fit is ",
      as.character(MSE_seq_new[optimal_ind], ".")
    )
  )
  print(plt)
}

calc_SSD <- function(X, r, f) {
  SSD_val <- purrr::map(1:nrow(X), ~ dist_euclidean(X[.x, ], f(r[.x, ]))^2) %>%
    unlist() %>%
    sum()
  SSD_val
}

update_parameterization <- function(
  time_points,
  r,
  X,
  f,
  n_knots,
  d,
  d2,
  gamma
) {
  full_r <- tidyr::expand_grid(time_points, r) %>%
    as.matrix(ncol = d + 1)
  new_param <- purrr::map(
    1:nrow(X),
    ~ projection_lpme(X[.x, ], f, full_r[.x, ], n_knots, d2, gamma)
  ) %>%
    purrr::reduce(cbind) %>%
    t()
  X_update <- purrr::map(1:nrow(full_r), ~ f(full_r[.x, ])) %>%
    purrr::reduce(rbind)

  list(
    parameterization = new_param,
    embedding = X_update
  )
}

calc_mse_cv <- function(
  leave_one_out,
  k,
  f,
  df,
  init_param,
  time_points,
  r,
  r_initial,
  n_knots,
  d,
  d_new,
  d_new2,
  D_out,
  D_coef,
  lambda,
  gamma,
  gamma2,
  r_full2,
  w,
  smoothing_method
) {
  if (leave_one_out == TRUE) {
    k <- length(time_points)
    folds <- sample(
      1:k,
      k,
      replace = FALSE
    )
  } else {
    folds <- sample(1:k, length(time_points), replace = TRUE)
  }
  cv_mse <- vector()
  r_mat <- gen_parameterization(r, n_knots, d)
  for (fold_idx in 1:k) {
    fold_times <- time_points[folds != fold_idx]
    r_full_cv <- tidyr::expand_grid(fold_times, r_mat) %>%
      as.matrix()

    x_vals <- purrr::map(1:nrow(r_full_cv), ~ f(r_full_cv[.x, ])) %>%
      purrr::reduce(rbind)
    spline_coefs <- list()
    for (idx in 1:length(fold_times)) {
      R <- cbind(rep(1, n_knots), r_mat)
      E <- calcE(r_mat, gamma)
      spline_coefs[[idx]] <- solve_spline(
        E,
        R,
        x_vals[x_vals[, 1] == fold_times[idx], -1],
        lambda[folds != fold_idx][idx],
        d,
        D_out
      ) %>%
        t() %>%
        matrix(nrow = 1)
    }
    coef_cv <- purrr::reduce(spline_coefs, rbind)

    if (smoothing_method == "spline") {
      r_cv <- matrix(
        r_full2[r_full2[, 1] %in% fold_times, ],
        ncol = 1
      )

      temp_param_cv <- cbind(rep(1, nrow(r_full2)), r_full2)
      E_cv <- calcE(r_cv, gamma2)
      param_cv <- temp_param_cv[temp_param_cv[, 2] %in% fold_times, ]

      sol_coef_cv <- solve_spline(
        E_cv,
        param_cv,
        coef_cv,
        w,
        1,
        D_coef
      )

      f_coef_cv <- function(t) {
        return_vec <- as.vector(
          t(sol_coef_cv[1:nrow(coef_cv), ]) %*%
            etaFunc(t, r_cv, gamma2) +
            t(sol_coef_cv[(nrow(coef_cv) + 1):(nrow(coef_cv) + 1 + 1), ]) %*%
              matrix(c(1, t), ncol = 1)
        )
        return_vec
      }

      f_new_cv <- function(t) {
        coefs <- f_coef_cv(t[1])
        coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
        return_vec <- t(coef_mat[1:n_knots, ]) %*%
          etaFunc(t[-1], r_initial, gamma) +
          t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*%
            matrix(c(1, t[-1]), ncol = 1)
        return(c(t[1], return_vec))
      }
    } else if (smoothing_method == "gp") {
      # invisible(
      #   capture.output({
      #     gp <- GPFDA::gpr(
      #       response = coef_cv,
      #       input = fold_times,
      #       Cov = "matern",
      #       meanModel = "t",
      #       nu = w
      #     )
      #   })
      # )
      #
      # f_new_cv <- function(t) {
      #   coefs <- GPFDA::gprPredict(train = gp, inputNew = t[1], noiseFreePred = TRUE)$pred.mean %>%
      #     as.vector()
      #   coef_mat <- matrix(coefs, n_knots + d + 1, byrow = TRUE)
      #   return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], r_initial, 4 - d) +
      #     t(coef_mat[(n_knots + 1):(n_knots + d + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
      #   c(t[1], return_vec)
      # }
    }

    temp_df <- df[!(df[, 1] %in% fold_times), ]
    temp_init_param <- init_param[!(df[, 1] %in% fold_times), ]
    cv_projections <- purrr::map(
      1:nrow(temp_df),
      ~ projection_lpme(temp_df[.x, ], f_new_cv, temp_init_param[.x, ])
    ) %>%
      purrr::reduce(rbind)

    cv_points <- purrr::map(
      1:nrow(cv_projections),
      ~ f_new_cv(cv_projections[.x, ])
    ) %>%
      purrr::reduce(cbind) %>%
      t()
    error_cv <- purrr::map(
      1:nrow(cv_points),
      ~ dist_euclidean(temp_df[.x, ], cv_points[.x, ])
    ) %>%
      purrr::reduce(c)
    cv_mse[fold_idx] <- mean(error_cv^2)
  }
  cv_mse
}

gen_parameterization <- function(r, n_knots, d) {
  r_bounds <- Rfast::colMinsMaxs(r)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = n_knots^(1 / d)
    )
  }
  r_mat <- as.matrix(expand.grid(r_list))
  r_mat
}

print_mse <- function(tuning, MSE) {
  print(
    paste0(
      "When gamma = ",
      as.character(tuning),
      ", MSD = ",
      as.character(MSE),
      "."
    )
  )
}
