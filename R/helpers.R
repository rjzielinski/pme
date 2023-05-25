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
