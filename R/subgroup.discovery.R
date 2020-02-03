.onUnload <- function (libpath) {
  library.dynam.unload("subgroup.discovery", libpath)
}

#' @useDynLib subgroup.discovery, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL
