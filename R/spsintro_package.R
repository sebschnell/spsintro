#' spsintro
#'
#' Code snippets and functions used for 'Introduction to Scientific Programming
#' and Simulation' course
#'
#' @name spsintro
#' @docType package
#' @importFrom Rcpp sourceCpp
#' @useDynLib spsintro

.onUnload <- function (libpath) {
  library.dynam.unload("spsintro", libpath)
}

