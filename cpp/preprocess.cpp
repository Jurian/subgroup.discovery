
#include <string>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
#include "preprocess.hpp"
#include "colworkers.hpp"

std::vector<arma::uword> sort_index(const Rcpp::NumericVector& v) {
  std::vector<arma::uword> index(v.size(), 0);
  for (arma::uword i = 0 ; i != index.size() ; i++) {
    index[i] = i;
  }
  std::sort(index.begin(), index.end(),
       [&](const arma::uword& a, const arma::uword& b) {
         return (v[a] < v[b]);
       }
  );
  return index;
}

std::map<arma::uword, Category> extractCategories(const Rcpp::IntegerVector& col) {
  std::map<arma::uword, Category> categories;

  const Rcpp::CharacterVector temp = col.attr("levels");
  std::string* levels = new std::string[temp.size()];
  for(arma::uword i = 0; i < temp.size(); i++) {
    levels[i] = temp[i];
  }

  for(arma::uword i = 0; i < col.size(); i++) {
    // Explicitly make sure we use an integer value here
    const arma::uword catID = col[i] - 1;
    if(categories.count(catID) == 0) {

      bool* categoryIndex = new bool[col.size()] {};
      categoryIndex[i] = true;

      Category category = {
        catID,
        levels[catID],
              categoryIndex
      };
      categories[catID] = category;
    } else {
      (categories.find(catID) -> second).index[i] = true;
    }
  }

  return categories;
}

PrimContainer preprocess(const Rcpp::DataFrame& df, const arma::vec& y, const double& alpha, const double& minSup, const bool& verbose) {

  const Rcpp::CharacterVector temp = df.names();

  std::string* colNames = new std::string[df.size()];
  ColumnType* colTypes = new ColumnType[df.size()];
  std::map<arma::uword, arma::uvec> orderMap;
  std::map<arma::uword, std::map<arma::uword, Category>> catMap;
  arma::mat M(df.nrows(), df.size());

  for(arma::uword col = 0; col < df.size(); col++) {
    colNames[col] = temp[col];

    switch(TYPEOF(df[col])){
    case REALSXP: {
      const Rcpp::NumericVector nv = df[col];
      for(arma::uword row = 0; row < nv.size(); row++) M(row, col) = nv[row];
      orderMap[col] = arma::uvec(sort_index(nv));
      colTypes[col] = NUMERIC;
      break;
    }
    case INTSXP: {

      if(Rf_isFactor(df[col])) {
        const Rcpp::IntegerVector nv = df[col];
        catMap[col] = extractCategories(nv);
        colTypes[col] = FACTOR;
        for(arma::uword row = 0; row < nv.size(); row++) M(row, col) = nv[row]-1;
      } else {
        const Rcpp::NumericVector nv = df[col];
        orderMap[col] = arma::uvec(sort_index(nv));
        colTypes[col] = NUMERIC;
        for(arma::uword row = 0; row < nv.size(); row++) M(row, col) = nv[row];
      }
      break;
    }
    default: {
      Rcpp::stop("Colmn %s is not numeric or factor.", colNames[col]);
    }
    }
  }

  PrimContainer pc = {
    colNames,
    colTypes,
    orderMap,
    catMap,
    M,
    y,
    alpha,
    minSup,
    verbose
  };

  return pc;
}
