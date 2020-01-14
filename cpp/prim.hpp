#ifndef _PRIM_HPP
#define _PRIM_HPP

// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
#include "preprocess.hpp"

arma::uvec peel(const PrimContainer& pc);

// [[Rcpp::export]]
arma::uvec prim(const Rcpp::DataFrame& df, const arma::vec& y, const double& quantile, const double& minSup, const bool& verbose = false);


#endif
