new_lpme <- function(msd,
                     sol_coef,
                     times,
                     r_init,
                     r_fit,
                     d,
                     D,
                     n_knots,
                     lambda,
                     gamma,
                     sol_coef_list,
                     r_fit_list) {
  lpme_list <- list(
    d = d,
    D = D,
    lambda = lambda,
    gamma = gamma,
    times = times,
    n_knots = n_knots,
    r_init = r_init,
    r_fit = r_fit,
    sol_coef = sol_coef,
    msd = msd,
    sol_coef_list = sol_coef_list,
    r_fit_list = r_fit_list
  )
  vctrs::new_vctr(lpme_list, class = "lpme")
}

is_lpme <- function(x) {
  inherits(x, "lpme")
}

lpme <- function(df,
                 d,
                 tuning_para_seq = c(0, exp(-15:10)),
                 alpha = 0.05,
                 max_comp = 500,
                 epsilon = 0.05,
                 max_iter = 100,
                 verbose = TRUE,
                 print_plots = TRUE,
                 SSD_ratio_threshold = 100,
                 increase_threshold = 1.05,
                 init = "full") {

  # Declare initial variable values ---------------------------------------
  time_points <- df[, 1] %>%
    unique()

  initialization <- initialize_lpme(df, init, time_points, d, alpha, max_comp)
  init_pme_list <- fit_init_pmes(df, time_points, init, initialization, d)
  spline_coef_list <- merge_spline_coefs(init_pme_list, d, time_points)

  coef_full <- spline_coef_list$coef_full
  x_test <- spline_coef_list$x_test
  r_full <- spline_coef_list$r_full
  r_full2 <- spline_coef_list$r_full2
  n_knots <- spline_coef_list$n_knots
  lambda <- spline_coef_list$lambda

  D_coef <- dim(coef_full)[2]
  D_new <- dim(x_test)[2] - 1
  D_new2 <- dim(x_test)[2]
  d_new <- dim(r_full)[2]
  d_new2 <- dim(r_full2)[2]
  n <- dim(x_test)[1]
  gamma <- 4 - d_new
  gamma2 <- 4 - d_new2
  X_new <- x_test
  I_new <- n
  t_initial <- r_full %>%
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

  for (tuning_ind in 1:length(tuning_para_seq)) {
    f_coef_list <- compute_f_coef(
      tuning_para_seq[tuning_ind],
      diag(weights),
      t_initial,
      r_full2,
      coef_full,
      d_new2,
      D_coef,
      gamma,
      gamma2
    )
    f_coef <- f_coef_list$f
    sol_coef <- f_coef_list$sol

    f_new <- function(t) {
      coefs <- f_coef(t[1])
      coef_mat <- matrix(coefs, n_knots + d_new + 1, byrow = TRUE)
      return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, gamma) +
        t(coef_mat[(n_knots + 1):(n_knots + d_new + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
      c(t[1], return_vec)
    }

    f0_new <- f_new

    full_t <- tidyr::expand_grid(time_points, t_initial) %>%
      as.matrix(ncol = d + 1)

    t_new2 <- purrr::map(
      1:nrow(X_new),
      ~projection_lpme(X_new[.x, ], f0_new, full_t[.x, ], n_knots, d_new, gamma)
    ) %>%
      purrr::reduce(cbind) %>%
      t()

    # X_initial_guess <- cbind(X_new, full_t)
    # projection.index.f0 <- function(x.init) {
    #   projection_lpme(
    #     x.init[1:D_new2],
    #     f0_new,
    #     x.init[(D_new2 + 1):(D_new2 + d_new + 1)],
    #     n_knots,
    #     d_new,
    #     gamma
    #   )
    # }

    # t_new2 <- matrix(
    #   t(apply(X_initial_guess, 1, projection.index.f0)),
    #   nrow = I_new
    # )

    x_fun <- purrr::map(
      1:I_new,
      ~ f_new(full_t[.x, ])
    ) %>%
      purrr::reduce(rbind)

    X_projection_index <- cbind(x_fun, full_t)

    SSD_new <- purrr::map(
      1:I_new,
      ~ dist_euclidean(
        X_new[.x, ],
        f_new(t_new2[.x, ])
      ) ^ 2
    ) %>%
      unlist() %>%
      sum()

    if (print_plots == TRUE) {
      plot_lpme(df, f_new, t_initial, d_new, D_new, time_points)
    }

    count <- 1
    SSD_ratio <- 10 * epsilon

    nearest_x <- calc_nearest_x(df, x_test)
    init_param <- calc_init_param(df, t_new2, nearest_x)

    data_initial <- cbind(df, init_param)

    cv_mse <- vector()
    for (time_idx in 1:length(time_points)) {
      r_bounds <- Rfast::colMinsMaxs(t_new2[, -1])
      r_list <- list()
      for (idx in 1:dim(r_bounds)[2]) {
        r_list[[idx]] <- seq(
          r_bounds[1, idx],
          r_bounds[2, idx],
          length.out = n_knots ^ (1 / d)
        )
      }
      # r_list[[1]] <- time_points
      r_mat <- as.matrix(expand.grid(r_list))
      r_full_cv <- tidyr::expand_grid(time_points[-time_idx], r_mat) %>%
        as.matrix()

      spline_coefs <- list()

      x_vals <- purrr::map(
        1:nrow(r_full_cv),
        ~ f_new(r_full_cv[.x, ])
      ) %>%
        purrr::reduce(rbind)

      for (idx in 1:length(time_points)) {
        if (time_idx == idx) {
          next
        }

        R <- cbind(rep(1, n_knots), r_mat)
        E <- calcE(r_mat, gamma)

        spline_coefs[[idx]] <- solve_spline(
          E,
          R,
          x_vals[x_vals[, 1] == time_points[idx], -1],
          lambda[idx],
          d_new,
          D_new
        ) %>%
          t() %>%
          matrix(nrow = 1)
      }

      coef_cv <- purrr::reduce(spline_coefs, rbind)
      x_test_cv <- x_vals
      X_new_cv <- x_test_cv
      # I_new <- length(theta_hat_new)
      I_new_cv <- nrow(X_new_cv)

      r_cv <- matrix(
        r_full2[r_full2[, 1] != time_points[time_idx], ],
        ncol = 1
      )

      T_new2 <- cbind(rep(1, nrow(r_full2)), r_full2)
      E_new_cv <- calcE(r_cv, gamma2)
      T_new_cv <- T_new2[T_new2[, 2] != time_points[time_idx], ]

      sol_coef_cv <- solve_spline(
        E_new_cv,
        T_new_cv,
        coef_cv,
        tuning_para_seq[tuning_ind],
        d_new2,
        D_coef
      )

      f_coef_cv <- function(t) {
        return_vec <- as.vector(
          # t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, t_new, gamma) +
          t(sol_coef_cv[1:nrow(coef_cv), ]) %*% etaFunc(t, r_cv, gamma2) +
            t(sol_coef[(nrow(coef_cv) + 1):(nrow(coef_cv) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
        )
        return(return_vec)
      }

      f_new_cv <- function(t) {
        coefs <- f_coef_cv(t[1])
        coef_mat <- matrix(coefs, n_knots + d_new + 1, byrow = TRUE)
        return_vec <- t(coef_mat[1:n_knots, ]) %*% etaFunc(t[-1], t_initial, gamma) +
          t(coef_mat[(n_knots + 1):(n_knots + d_new + 1), ]) %*% matrix(c(1, t[-1]), ncol = 1)
        return(c(t[1], return_vec))
      }

      temp_data_initial <- data_initial[data_initial[, 1] == time_points[time_idx], ]
      proj_para_cv <- purrr::map(
        1:nrow(temp_data_initial),
        ~ try({
          projection_lpme(
            temp_data_initial[.x, 1:D_new2],
            f_new_cv,
            temp_data_initial[.x, (D_new2 + 1):(D_new2 + d_new + 1)]
          ) %>%
            t()
        })
      )

      errors <- vector()
      for (i in 1:length(proj_para_cv)) {
        errors[i] <- sum(class(proj_para_cv[[i]]) == "try-error") > 0
      }

      if (sum(errors) > 0) {
        proj_para_cv <- purrr::reduce(proj_para_cv, rbind)[-which(errors), ]
      } else {
        proj_para_cv <- purrr::reduce(proj_para_cv, rbind)
      }

      proj_points_cv <- t(apply(proj_para_cv, 1, f_new_cv))

      if (sum(errors) > 0) {
        temp_data_initial <- temp_data_initial[-which(errors), ]
      }

      proj_error_cv <- purrr::map(
        1:nrow(proj_points_cv),
        ~dist_euclidean(temp_data_initial[.x, 1:D_new2], proj_points_cv[.x, ])
      ) %>%
        purrr::reduce(c)
      # proj_error_cv <- dist_euclideanC_vec(
      #   temp_data_initial[, 1:D_new2],
      #   proj_points_cv
      # ) %>%
      #   as.vector()
      cv_mse[time_idx] <- mean(proj_error_cv ^ 2)
    }

    df_n <- sapply(time_points, function(x) nrow(df[df[, 1] == x, ]))
    MSE_new <- stats::weighted.mean(cv_mse, df_n)
    MSE_seq_new[tuning_ind] <- MSE_new

    if (verbose == TRUE) {
      print(
        paste(
          "When gamma = ",
          as.character(tuning_para_seq[tuning_ind]),
          ", MSD = ",
          as.character(MSE_new),
          "."
        )
      )
    }

    SOL_coef[[tuning_ind]] <- sol_coef
    TNEW_new[[tuning_ind]] <- t_new2
    coefs[[tuning_ind]] <- r_mat
    x_funs[[tuning_ind]] <- x_fun
    functions[[tuning_ind]] <- f_new
    func_coef[[tuning_ind]] <- f_coef
  }

  optimal_ind <- min(which(MSE_seq_new == min(MSE_seq_new)))
  sol_opt <- SOL_coef[[optimal_ind]]
  t_new_opt <- TNEW_new[[optimal_ind]]
  coefs_opt <- coefs[[optimal_ind]]
  f.optimal <- functions[[optimal_ind]]

  if (verbose == TRUE) {
    plot(
      log(tuning_para_seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      xlab = "Log Gamma",
      ylab = "MSD",
      type = "l"
    )
    graphics::lines(
      log(tuning_para_seq[1:length(MSE_seq_new)]),
      MSE_seq_new,
      type = "p",
      pch = 20,
      col = "orange",
      cex = 2
    )
    graphics::abline(
      v = log(tuning_para_seq[optimal_ind]),
      lwd = 1.5,
      col = "darkgreen",
      lty = 2
    )
    print(
      paste(
        "The optimal tuning parameter is ",
        as.character(tuning_para_seq[optimal_ind]),
        ", and the MSD of the optimal fit is ",
        as.character(MSE_seq_new[optimal_ind], ".")
      )
    )
  }

  lpme_out <- new_lpme(
    msd = MSE_seq_new,
    sol_coef = sol_opt,
    times = as.vector(time_points),
    r_init = t_initial,
    r_fit = t_new_opt,
    d = d_new,
    D = D_new2,
    n_knots = n_knots,
    lambda = lambda,
    gamma = gamma,
    sol_coef_list = SOL_coef,
    r_fit_list = TNEW_new
  )

  lpme_out
}

embed <- function(object, x) {
  n_times <- length(object$times)
  nu <- 4 - object$d
  f_coef <- function(t) {
    vec <- as.vector(
      t(object$sol_coef[1:n_times, ]) %*% etaFunc(t, as.matrix(object$times, ncol = 1), 3) +
        t(object$sol_coef[(n_times + 1):(n_times + 2), ]) %*% matrix(c(1, t), ncol = 1)
    )
    return(vec)
  }

  coefs <- f_coef(x[1])
  coef_mat <- matrix(coefs, object$n_knots + object$d + 1, byrow = TRUE)
  vec <- as.vector(
    t(coef_mat[1:object$n_knots, ]) %*% etaFunc(x[-1], object$r_init, nu) +
      t(coef_mat[(object$n_knots + 1):(object$n_knots + object$d + 1), ]) %*% matrix(c(1, x[-1]), ncol = 1)
  )
  return(c(x[1], vec))
}

parameterize <- function(object, x) {
  embeddings <- purrr::map(
    1:nrow(object$r_fit),
    ~ embed(object, object$r_fit[.x, ])
  ) %>%
    purrr::reduce(rbind)
  nearest <- calc_nearest_x(matrix(x, nrow = 1), embeddings[embeddings[, 1] == x[1], ])
  estimate <- projection_lpme(x, function(t) embed(object, t), object$r_fit[object$r_fit[, 1] == x[1], ][nearest, ])

  return(estimate)
}

initialize_lpme <- function(df, init, time_points, d, alpha, max_comp) {
  init_timevals <- list()
  init_theta_hat <- list()
  init_centers <- list()
  init_sigma <- list()
  init_clusters <- list()

  if (init %in% c("first", "full")) {
    if (init == "first") {
      init_df <- df[df[, 1] == time_points[1], -1]
    } else if (init == "full") {
      init_dimension_size <- dim(df[, -1])
      init_D <- init_dimension_size[2]
      init_n <- init_dimension_size[1]
      lambda <- 4 - d
      init_N0 <- 20 * init_D
      for (idx in 1:length(time_points)) {
        init_df_temp <- df[df[, 1] == time_points[idx], -1]
        init_est_temp <- hdmde(init_df_temp, init_N0, alpha, max_comp)
        init_timevals[[idx]] <- rep(time_points[idx], dim(init_est_temp$mu)[1])
        init_theta_hat[[idx]] <- init_est_temp$theta_hat
        init_centers[[idx]] <- init_est_temp$mu
        init_sigma[[idx]] <- init_est_temp$sigma
        init_clusters[[idx]] <- init_est_temp$km
      }
    }

    init_timevals <- purrr::reduce(init_timevals, c)
    init_theta_hat <- purrr::reduce(init_theta_hat, c)
    init_centers <- purrr::reduce(init_centers, rbind)
    init_sigma <- purrr::reduce(init_sigma, c)
    init_W <- diag(init_theta_hat)
    init_X <- init_centers
    init_I <- length(init_theta_hat)

    init_dissimilarity_matrix <- as.matrix(stats::dist(init_X))
    init_isomap <- vegan::isomap(init_dissimilarity_matrix, ndim = d, k = 10)

    init_list <- list(
      timevals = init_timevals,
      theta_hat = init_theta_hat,
      centers = init_centers,
      sigma = init_sigma,
      clusters = init_clusters,
      W = init_W,
      X = init_X,
      I = init_I,
      param = init_isomap
    )
    return(init_list)
  } else {
    print("Please choose `init` option 'first' or 'full'.")
  }
}

fit_init_pmes <- function(df, time_points, init, initialization, d) {
  pme_results <- list()
  pme_coefs <- list()
  pme_coefs2 <- list()
  funcs <- list()
  clusters <- list()
  centers <- list()
  x_test <- list()
  r <- list()
  r2 <- list()
  num_clusters <- rep(0, length(time_points))
  errors <- vector()

  init_timevals <- initialization$timevals
  init_theta_hat <- initialization$theta_hat
  init_clusters <- initialization$clusters
  init_centers <- initialization$centers
  init_sigma <- initialization$sigma
  init_isomap <- initialization$param

  for (idx in 1:length(time_points)) {
    df_temp <- df[df[, 1] == time_points[idx], ]
    if (init %in% c("first", "full")) {
      pme_results[[idx]] <- pme(
        x_obs = df_temp[, -1],
        d = d,
        tuning_para_seq = exp(-20:10),
        initialization = list(
          parameterization = matrix(
            init_isomap$points[init_timevals == time_points[idx], ],
            nrow = length(init_theta_hat[init_timevals == time_points[idx]])
          ),
          theta_hat = init_theta_hat[init_timevals == time_points[idx]],
          centers = init_centers[init_timevals == time_points[idx], ],
          sigma = init_sigma[idx],
          km = init_clusters[[idx]]
        ),
        verbose = FALSE
      )
    } else {
      pme_results[[idx]] <- pme(
        x_obs = df_temp[, -1],
        d = d,
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
    pme_coefs[[idx]] <- pme_results[[idx]]$coefs[[opt_run]][1:num_clusters[idx], ]
    pme_coefs2[[idx]] <- pme_results[[idx]]$coefs[[opt_run]][(num_clusters[idx] + 1):(num_clusters[idx] + d + 2), ] %>%
      t() %>%
      matrix(nrow = 1)
    x_test[[idx]] <- apply(
      pme_results[[idx]]$params[[opt_run]],
      1,
      funcs[[idx]]
    ) %>%
      t()
    r[[idx]] <- pme_results[[idx]]$params[[opt_run]]
    r2[[idx]] <- time_points[idx]
    errors[idx] <- pme_results[[idx]]$MSD[opt_run]
  }

  init_pme_list <- list(
    pme_results = pme_results,
    pme_coefs = pme_coefs,
    pme_coefs2 = pme_coefs2,
    funcs = funcs,
    clusters = clusters,
    centers = centers,
    x_test = x_test,
    r = r,
    r2 = r2,
    num_clusters = num_clusters,
    errors = errors
  )
  init_pme_list
}

merge_spline_coefs <- function(pme_list, d, time_points) {
  lambda <- vector()
  x_vals <- list()
  spline_coefs <- list()
  gamma <- 4 - d

  r_full <- purrr::reduce(pme_list$r, rbind)
  length_r <- ceiling(max(pme_list$num_clusters)^(1 / d))
  r_bounds <- Rfast::colMinsMaxs(r_full)
  r_list <- list()
  for (idx in 1:dim(r_bounds)[2]) {
    r_list[[idx]] <- seq(
      r_bounds[1, idx],
      r_bounds[2, idx],
      length.out = length_r
    )
  }
  r_full <- as.matrix(expand.grid(r_list))
  n_knots <- nrow(r_full)
  r_full2 <- purrr::reduce(pme_list$r2, rbind)

  for (time_idx in 1:length(time_points)) {
    lambda[time_idx] <- pme_list$pme_results[[time_idx]]$tuning
    f <- pme_list$funcs[[time_idx]]
    output <- purrr::map(
      1:n_knots,
      ~f(r_full[.x, ])
    ) %>%
      purrr::reduce(rbind)
    x_vals[[time_idx]] <- cbind(time_points[time_idx], output)

    R <- cbind(rep(1, n_knots), r_full)
    E <- calcE(r_full, gamma)

    spline_coefs[[time_idx]] <- solve_spline(E, R, output, lambda[time_idx], ncol(r_full), ncol(output)) %>%
      t() %>%
      matrix(nrow = 1)
  }
  coef_full <- purrr::reduce(spline_coefs, rbind)
  x_test <- purrr::reduce(x_vals, rbind)

  out_list <- list(
    r_full = r_full,
    r_full2 = r_full2,
    coef_full = coef_full,
    x_test = x_test,
    n_knots = n_knots,
    lambda = lambda
  )
  out_list
}

compute_f_coef <- function(w_new, W, t_new, r_full2, coef_full, d_new2, D_coef, gamma, gamma2) {
  T_new <- cbind(rep(1, nrow(t_new)), t_new)
  T_new2 <- cbind(rep(1, nrow(r_full2)), r_full2)

  E_new <- calcE(t_new, gamma)
  E_new2 <- calcE(r_full2, gamma2)

  sol_coef <- solve_weighted_spline(E_new2, W, T_new2, coef_full, w_new, d_new2, D_coef)
  f_coef <- function(t) {
    return_vec <- as.vector(
      t(sol_coef[1:nrow(coef_full), ]) %*% etaFunc(t, r_full2, gamma2) +
        t(sol_coef[(nrow(coef_full) + 1):(nrow(coef_full) + d_new2 + 1), ]) %*% matrix(c(1, t), ncol = 1)
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
    ~f(pred_grid[.x, ])
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

