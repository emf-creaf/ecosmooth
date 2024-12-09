#' ecosmooth: Kernel-based smoothing of multivariate ecological data
#'
#' Performs smoothing of multivariate ecological data on the basis of Gaussian kernels
#'
#' @name ecosmooth-package
#' @aliases ecosmooth ecosmooth-package
#' @docType package
#' @author \strong{Maintainer}: Miquel De Cáceres
#' \email{miquelcaceres@@gmail.com}
#' [ORCID](https://orcid.org/0000-0001-7132-2080)
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.dist dist
#' @useDynLib ecosmooth, .registration = TRUE
## usethis namespace: end
NULL
