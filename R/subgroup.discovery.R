# Packages depending on RcppParallel often have issues with ASAN / UBSANwarnings,
# originating from use of RcppParallel's bundled TBB library
# To alleviate this, with RcppParallel 5.0.0, it is now possible to configure the
# parallel backend used by RcppParallel by setting the environment variable.
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

.onUnload <- function (libpath) {
  library.dynam.unload("subgroup.discovery", libpath)
}

#' @useDynLib subgroup.discovery, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
NULL
