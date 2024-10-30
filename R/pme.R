#' Principal Manifold Estimation
#'
#' This function still requires completed documentation
#'
#' @param data A numeric matrix of high-dimensional data.
#' @param d A positive integer representing the intrinsic dimension.
#' @param initialization A list of values providing an initialization for `pme()`. It is not recommended to supply these values manually.
#' @param initialization_algorithm The name of a manifold learning algorithm to use for finding the initial parameterizations. Options include "isomap", "diffusion_maps", and "laplacian_eigenmaps".
#' @param initialization_type Choose whether to use cluster centers or represent clusters by subsampling from the associated points. Options: "centers" or "subsample".
#' @param lambda A vector of smoothing values to be considered.
#' @param alpha The significant level to be used when testing the need for additional clusters in data reduction.
#' @param min_clusters The minimum number of clusters allowed in data reduction.
#' @param max_clusters The maximum number of clusters allowed in data reduction.
#' @param epsilon Threshold for change in SSD to stop iterations for a given smoothing value.
#' @param max_iter The maximum number of iterations allowed for a given smoothing value.
#' @param SSD_ratio_threshold  The maximum increase in SSD allowed before stopping iterations for a given smoothing value.
#' @param print_plots A logical value indicating whether plots of the estimated manifolds should be produced.
#' @param verbose A logical value indicating whether updates should be printed.
#'
#' @return An object of type pme.
#' @export
pme <- function(data,
    d,
    initialization = NULL,
    initialization_algorithm = "isomap",
    initialization_type = "centers",
    lambda = exp(-15:5),
    alpha = 0.01,
    min_clusters = 0,
    max_clusters = 100,
    epsilon = 0.05,
    max_iter = 100,
    SSD_ratio_threshold = 5,
    print_plots = FALSE,
    verbose = FALSE) {

  # Initial variable assignments --------------------------
  n <- dim(data)[1]
  D <- dim(data)[2]

  if (min_clusters == 0) {
    min_clusters <- 5 * d
  }
  if (min_clusters > max_clusters) {
    max_clusters <- nrow(data) - 1
  }

  # Initialization ----------------------------------------
  if (is.null(initialization)) {
    if (initialization_type == "subsample") {
      initialization <- initialize_pme(
        data,
        d,
        min_clusters,
        alpha,
        max_clusters,
        algorithm = initialization_algorithm,
        component_type = "subsample"
      )
    } else {
      initialization <- initialize_pme(
        data,
        d,
        min_clusters,
        alpha,
        max_clusters,
        algorithm = initialization_algorithm
      )
    }
  }
  weights <- diag(initialization$theta_hat)
  X <- initialization$centers
  I <- length(initialization$theta_hat)

  mse <- vector()
  coefs <- list()
  parameterization <- list()
  embeddings <- list()

  for (tuning_idx in 1:length(lambda)) {
    spline_coefs <- calc_coefficients(
      X,
      initialization$parameterization,
      weights,
      lambda[tuning_idx]
    )

    params <- initialization$parameterization

    f_embedding <- function(parameters) {
      as.vector(
        (t(spline_coefs[1:I, ]) %*% etaFunc(parameters, params, 4 - d)) +
          (t(spline_coefs[(I + 1):(I + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
      )
    }

    f0 <- f_embedding

    params <- calc_params(f_embedding, X, params)

    SSD <- calc_SSD(f_embedding, X, params)

    count <- 1
    SSD_ratio <- 10 * epsilon

    if (print_plots == TRUE) {
      plot_pme(f_embedding, data, X, spline_coefs, params, d)
    }

    while ((SSD_ratio > epsilon) & (SSD_ratio <= SSD_ratio_threshold) & (count <= (max_iter - 1))) {
      SSD_prev <- SSD
      f0 <- f_embedding
      params_prev <- params
      coefs_prev <- spline_coefs

      spline_coefs <- calc_coefficients(
        X,
        params,
        weights,
        lambda[tuning_idx]
      )

      f_embedding <- function(parameters) {
        as.vector(
          (t(spline_coefs[1:I, ]) %*% etaFunc(parameters, params, 4 - d)) +
            (t(spline_coefs[(I + 1):(I + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
        )
      }

      params <- calc_params(f_embedding, X, params)

      SSD <- calc_SSD(f_embedding, X, params)

      SSD_ratio <- abs(SSD - SSD_prev) / SSD_prev
      count <- count + 1

      if (print_plots == TRUE) {
        plot_pme(f_embedding, data, X, spline_coefs, params, d)
      }


      if (SSD_ratio > SSD_ratio_threshold) {
        f_embedding <- f0
        params <- params_prev
        spline_coefs <- coefs_prev
        SSD <- SSD_prev
      }

      if (verbose == TRUE) {
        print_SSD(lambda[tuning_idx], SSD, SSD_ratio, count)
      }
    }

    if (print_plots == TRUE) {
      plot_pme(f_embedding, data, X, spline_coefs, params, d)
    }

    mse[tuning_idx] <- calc_msd(data, initialization$km, f_embedding, params, D, d)
    if (verbose == TRUE) {
      print(paste("When lambda = ", as.character(lambda[tuning_idx]), ", MSD = ", as.character(mse[tuning_idx]), "."))
    }
    coefs[[tuning_idx]] <- spline_coefs
    parameterization[[tuning_idx]] <- params
    embeddings[[tuning_idx]] <- function(parameters) {
      as.vector(
        (t(coefs[[tuning_idx]][1:I, ]) %*% etaFunc(parameters, parameterization[[tuning_idx]], 4 - d)) +
          (t(coefs[[tuning_idx]][(I + 1):(I + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
      )
    }

    if (tuning_idx >= 4) {
      if (!is.unsorted(mse[(tuning_idx - 3):tuning_idx])) {
        break
      }
    }
  }

  optimal_idx <- min(which(mse == min(mse)))
  coefs_opt <- coefs[[optimal_idx]]
  params_opt <- parameterization[[optimal_idx]]
  embedding_opt <- function(parameters) {
    as.vector(
      (t(coefs_opt[1:I, ]) %*% etaFunc(parameters, params_opt, 4 - d)) +
        (t(coefs_opt[(I + 1):(I + d + 1), ]) %*% matrix(c(1, parameters), ncol = 1))
    )
  }

  if (verbose == TRUE) {
    paste0(
      "The optimal tuning parameter is ",
      as.character(lambda[optimal_idx]),
      ", and the MSD of the optimal fit is ",
      as.character(mse[optimal_idx]),
      "."
    )
  }

  pme_out <- new_pme(
    # embedding_map = embeddings[[optimal_idx]],
    embedding_map = embedding_opt,
    knots = initialization$km,
    knot_weights = initialization$theta_hat,
    kernel_coefs = coefs_opt[1:I, ],
    polynomial_coefs = coefs_opt[(I + 1):(I + d + 1), ],
    tuning = lambda[optimal_idx],
    MSD = mse,
    coefs = coefs,
    parameterization = parameterization,
    tuning_vec = lambda,
    embeddings = embeddings,
    initialization_algorithm = initialization_algorithm
  )
  pme_out
}

#' Create New PME Object
#'
#' @param embedding_map A function mapping from low-dimensional to
#' high-dimensional space.
#' @param knots A numeric matrix of mixture components.
#' @param knot_weights A numeric vector of the weights associated with
#' each mixture component.
#' @param kernel_coefs A numeric matrix of the coefficients in the first
#' part of a spline function.
#' @param polynomial_coefs A numeric matrix of the coefficients in the second
#' part of a spline function.
#' @param tuning A numeric value describing the smoothing parameter.
#' @param MSD A numeric vector of the mean squared distances associated with
#' the estimated embedding maps for each smoothing value.
#' @param coefs A list of the estimated spline coefficients of the embedding
#' maps for each smoothing value.
#' @param parameterization A list of the estimated parameterizations of the mixture
#' components for each smoothing value.
#' @param tuning_vec A numeric vector of smoothing values.
#' @param embeddings A list of the embedding maps estimated at each smoothing
#' value.
#' @param initialization_algorithm A character value indicating the algorithm used to find the initial parameterization.
#'
#' @return An object of type PME.
#'
#' @noRd
new_pme <- function(embedding_map,
                    knots,
                    knot_weights,
                    kernel_coefs,
                    polynomial_coefs,
                    tuning,
                    MSD,
                    coefs,
                    parameterization,
                    tuning_vec,
                    embeddings,
                    initialization_algorithm) {
  pme_list <- list(
    embedding_map = embedding_map,
    knots = knots,
    knot_weights = knot_weights,
    kernel_coefs = kernel_coefs,
    polynomial_coefs = polynomial_coefs,
    tuning = tuning,
    MSD = MSD,
    coefs = coefs,
    parameterization = parameterization,
    tuning_vec = tuning_vec,
    embeddings = embeddings,
    initialization = initialization_algorithm
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

#' Create PME Initialization
#'
#' @param x A numeric matrix of data.
#' @param d The intrinsic dimension.
#' @param min_clusters The minimum number of mixture components.
#' @param alpha Significance level.
#' @param max_clusters Maximum number of components.
#'
#' @return A list used to initialize PME.
#'
#' @noRd
initialize_pme <- function(x, d, min_clusters, alpha, max_clusters, component_type = "centers", algorithm = "isomap", subsample_size = 5) {
  est <- hdmde_mod(x, min_clusters, alpha, max_clusters)
  if (component_type == "subsample") {
    cluster_points <- matrix(nrow = 1, ncol = ncol(x))
    point_weights <- vector()
    for (cluster in 1:nrow(est$mu)) {
      temp_x <- x[est$km$cluster == cluster, ] %>%
        matrix(ncol = ncol(x))
      # cluster_distances <- map(
      #   1:nrow(temp_x),
      #   ~ dist_euclidean(est$mu[cluster, ], temp_x[.x, ])
      # ) %>%
      #   reduce(c)
      cluster_sample <- sample(
        1:nrow(temp_x),
        size = subsample_size,
        replace = TRUE
      )
      points <- unique(cluster_sample)
      n_occurences <- table(cluster_sample)

      cluster_points <- rbind(
        cluster_points,
        temp_x[points, ]
      )
      point_weights <- c(
        point_weights,
        est$theta_hat[cluster] * n_occurences
      )
    }
    cluster_points <- cluster_points[-1, ]
    est_order <- order(cluster_points[, 1])
    centers <- cluster_points[est_order, ]
    W <- diag(point_weights)
    theta_hat <- point_weights
  } else {
    est_order <- order(est$mu[, 1])
    theta_hat <- est$theta_hat[est_order]
    centers <- est$mu[est_order, ]
    W <- diag(theta_hat)
  }
  sigma <- est$sigma
  X <- centers
  I <- nrow(centers)

  if (algorithm == "diffusion_maps") {
    init_parameterization <- dimRed::embed(
      X,
      "DiffusionMaps",
      ndim = d,
      .mute = c("message", "output")
    )
  # } else if (algorithm == "hessian_eigenmaps") {
  #   init_parameterization <- dimRed::embed(
  #     X,
  #     "HLLE",
  #     knn = floor(sqrt(nrow(X))),
  #     ndim = d
  #     # .mute = c("message", "output")
  #   )
  } else if (algorithm == "laplacian_eigenmaps") {
    init_parameterization <- dimRed::embed(
      X,
      "LaplacianEigenmaps",
      knn = floor(sqrt(nrow(X))),
      ndim = d,
      .mute = c("message", "output")
    )
  # } else if (algorithm == "lle") {
  #   init_parameterization <- dimRed::embed(
  #     X,
  #     "LLE",
  #     knn = floor(sqrt(nrow(X))),
  #     ndim = d
  #     # .mute = c("message", "output")
  #   )
  } else {
    init_parameterization <- dimRed::embed(
      X,
      "Isomap",
      knn = floor(sqrt(nrow(X))),
      ndim = d,
      .mute = c("message", "output")
    )
    # dissimilarity <- as.matrix(stats::dist(X))
    # init_parameterization <- vegan::isomap(dissimilarity, ndim = d, k = floor(sqrt(nrow(dissimilarity))))
  }

  output <- dimRed::as.data.frame(init_parameterization) %>%
    as.matrix()

  list(
    parameterization = matrix(output[, 1:d], nrow = nrow(X)),
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
#' @param weights Numeric matrix of cluster weights.
#' @param w smoothing parameter.
#'
#' @return A matrix of spline coefficients.
#'
#' @noRd
calc_coefficients <- function(X, t, weights, w) {
  t_val <- cbind(rep(1, nrow(t)), t)
  E <- calcE(t, 4 - ncol(t))
  solve_weighted_spline(E, weights, t_val, X, w, ncol(t), ncol(X))
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
calc_params <- function(f, X, init_params) {
  params <- purrr::map(1:nrow(X), ~ projection_pme(X[.x, ], f, init_params[.x, ])) %>%
    unlist() %>%
    matrix(nrow = nrow(X), byrow = TRUE)
 params
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
  SSD_val <- purrr::map(1:nrow(X), ~ dist_euclidean(X[.x, ], f(t[.x, ]))^2) %>%
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
print_SSD <- function(tuning_val, SSD_new, SSD_ratio, count) {
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

#' Plot a PME Object
#'
#' @param f A function of the embedding map.
#' @param x A numeric matrix of the unreduced data.
#' @param centers A numeric matrix of the mixture component centers.
#' @param sol A numeric matrix of the embedding map coefficients.
#' @param t A numeric matrix of the parameterization of the component centers.
#' @param d The intrinsic dimension.
#'
#' @noRd
plot_pme <- function(f, x, centers, sol, t, d) {
  # pred_grid <- calc_tnew(centers, t, sol, I, d, lambda)
  pred_grid <- calc_params(f, centers, t)
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
    ~ f(unlist(as.vector(pred_grid[.x, ]))) %>%
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
        opacity = 0.15,
        marker = list(size = 3)
      )
    print(plt)
  }
}

#' Calculate Mean Squared Distance
#'
#' @param x A numeric matrix of data
#' @param km A kmeans object describing reduced data.
#' @param f An embedding map.
#' @param t A numeric matrix of parameterizations for the mixture components
#' described in `km`.
#' @param D The dimension of the original dataset.
#' @param d The intrinsic dimension.
#'
#' @return A numeric value.
#'
#' @noRd
calc_msd <- function(x, km, f, t, D, d) {
  data_initial <- matrix(0, nrow = 1, ncol = D + d)
  center_order <- order(km$centers[, 1])
  for (i in 1:max(km$cluster)) {
    index_temp <- which(km$cluster == center_order[i])
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
    ~ projection_pme(data_initial[.x, 1:D], f, data_initial[.x, (D + 1):(D + d)]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)
  proj_points <- purrr::map(
    1:nrow(proj_para),
    ~ f(proj_para[.x, ]) %>%
      t()
  ) %>%
    purrr::reduce(rbind)

  mse <- purrr::map(
    1:nrow(data_initial),
    ~ dist_euclidean(data_initial[.x, 1:D], proj_points[.x, ])^2
  ) %>%
    unlist() %>%
    mean()

  mse
}
