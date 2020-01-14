#ifndef _PREPROCESS_HPP
#define _PREPROCESS_HPP

// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <map>
#include <RcppArmadillo.h>

enum ColumnType { NUMERIC, FACTOR };

struct Category {
  arma::uword factor;
  std::string level;
  bool* index;
};

struct PrimContainer {
  std::string* colNames;
  ColumnType* colTypes;
  std::map<arma::uword, arma::uvec> orderMap;
  std::map<arma::uword, std::map<arma::uword, Category>> catMap;
  arma::mat M;
  arma::vec y;
  double alpha;
  double minSup;
  bool verbose;
};

std::vector<arma::uword> sort_index(const Rcpp::NumericVector& v);

std::map<arma::uword, Category> extractCategories(const Rcpp::IntegerVector& col);

PrimContainer preprocess(const Rcpp::DataFrame& df, const arma::vec& y, const double& alpha, const double& minSup, const bool& verbose = false);


#endif
