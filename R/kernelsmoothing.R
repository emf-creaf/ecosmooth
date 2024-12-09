
#' Kernel-based smoothing of multivariate data
#'
#' @param x A rectangular matrix (with observations in rows and variables in columns) or an object of class \code{\link{dist}} (with resemblance values between observations).
#' @param dist_kernel An object of class \code{\link{dist}} containing distances in the kernel space.
#' @param kernel_scale Bandwidth of the Gaussian kernel.
#' @param add Boolean flag to set to zero negative squared distances in the distance-based smoothing (otherwise corresponding values will be NA).
#'
#' @return
#' \itemize{
#'   \item{Function \code{kernelsmoothing} returns either a rectangular matrix (with observations in rows and variables in columns) or an object of class \code{\link{dist}} (with resemblance values between observations),
#' depending on the input \code{x}.}
#'   \item{Function \code{kernelmatrix} returns a square numeric matrix with kernel values.}
#' }
#' @export
#' @name kernelsmoothing
#' @examples
#' # Create original data
#' x <- matrix(rnorm(100), nrow = 5)
#'
#' # Original distance matrix
#' d <- dist(x)
#'
#' # Kernel matrix based on d space (bandwidth = 2)
#' k <- kernelmatrix(d, 2)
#' k
#'
#' # Self-smoothing of distance data (bandwidth = 2)
#' kernelsmoothing(d, d, 2)
#'
#' # Check equivalence with self-smoothing of rectangular data
#' dist(kernelsmoothing(x,d,2))
kernelsmoothing <- function(x, dist_kernel, kernel_scale, add = TRUE) {
  if(!inherits(x, "matrix") && !inherits(x, "dist")) stop("x has to be a numeric matrix or an object of class dist")
  umat <- kernelmatrix(dist_kernel, kernel_scale, normalize = TRUE)

  if(inherits(x, "matrix")) {
    return(.rectangularMatrixSmoothing(x, umat))
  } else {
    return(as.dist(.distanceMatrixSmoothing(as.matrix(x), umat, add)))
  }
}

#' @rdname kernelsmoothing
#' @param normalize A boolean flag to normalize kernel values to sum one for each row
#' @export
kernelmatrix <- function(dist_kernel, kernel_scale, normalize = TRUE) {
  if(!inherits(dist_kernel, "dist")) stop("dist_kernel has to be an object of class dist")
  if(!is.numeric(kernel_scale))  stop("kernel_scale has to be a positive numeric value")
  if(kernel_scale <= 0)  stop("kernel_scale has to be a positive numeric value")
  dist_kernel <- as.matrix(dist_kernel)
  N <- nrow(dist_kernel)
  umat_kernel <- matrix(NA, N, N)
  for(i1 in 1:N) {
    for(i2 in 1:N) {
      umat_kernel[i1, i2] <- exp(-1*(dist_kernel[i1,i2]^2)/(2*(kernel_scale^2)))
    }
  }
  if(normalize) {
    for(i in 1:N) umat_kernel[i,] <- umat_kernel[i, ]/sum(umat_kernel[i, ])
  }
  return(umat_kernel)
}

