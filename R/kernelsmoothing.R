
#' Kernel-based smoothing of multivariate data
#'
#' @param x A rectangular matrix (with observations in rows and variables in columns) or an object of class \code{\link{dist}} (with resemblance values between observations).
#' @param dist_kernel An object of class \code{\link{dist}} containing distances in the kernel space.
#' @param kernel_scale Bandwidth of the Gaussian kernel.
#'
#' @return Either a rectangular matrix (with observations in rows and variables in columns) or an object of class \code{\link{dist}} (with resemblance values between observations),
#' depending on the input \code{x}.
#' @export
#'
#' @examples
kernelsmoothing <- function(x, dist_kernel, kernel_scale) {
  if(!inherits(x, ""))
  print("Hello, world!")
}
