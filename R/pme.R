
#' Principal Manifold Estimation
#'
#' This function still requires completed documentation
#'
#' @param x_obs Value
#' @param d Value
#' @param initialization Value
#' @param N0 Value
#' @param tuning.para.seq Value
#' @param alpha Value
#' @param max.comp Value
#' @param epsilon Value
#' @param max.iter Value
#' @param SSD_ratio_threshold Value
#' @param print_plots Value
#' @param verbose Value
#'
#' @return An object of type pme.
#' @export
#'
pme <- function(x_obs, d, initialization = NULL, N0 = 20 * D, tuning.para.seq = exp(-15:5), alpha = 0.05, max.comp = 100, epsilon = 0.05, max.iter = 100, SSD_ratio_threshold = 10, print_plots = FALSE, verbose = TRUE) {

  ### Start of pme_prelim()
  n <- dim(x_obs)[1]
  D <- dim(x_obs)[2]
  lambda <- 4 - d

  if (N0 == 0) {
    N0 <- 20 * D
  }
  if (N0 > max.comp) {
    max.comp <- 2 * N0
  }

  if (is.null(initialization)) {
    ### Start of initialize_pme()
    est <- hdmde(x_obs, N0, alpha, max.comp)
  }

}

new_pme <- function() {

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


