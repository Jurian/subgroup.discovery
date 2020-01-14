#ifndef _MAPREDUCE_HPP
#define _MAPREDUCE_HPP

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "preprocess.hpp"

enum BoxType { NUM_LEFT, NUM_RIGHT, CATEGORY };

struct SubBox {
  std::vector<arma::uword> remove;
  std::vector<arma::uword> keep;
  arma::uword col;
  std::string value;
  double quality;
  BoxType side;

  int compare(const SubBox& cmp) const;
};

struct ColWorker : public RcppParallel::Worker {

  const arma::mat M;
  const arma::vec y;
  const double p, minSup, N;
  const std::map<arma::uword, arma::uvec> orderMap;
  const std::map<arma::uword, std::map<arma::uword, Category>> catMap;
  const std::string* colNames;
  const ColumnType* colTypes;

  const arma::uword masked;
  const bool* mask;

  bool subBoxFound;
  SubBox bestSubBox;

  // constructors
  ColWorker (const PrimContainer& pc, const arma::uword masked, const bool* mask)
    : M(pc.M),
      y(pc.y),
      p(pc.alpha),
      minSup(pc.minSup),
      N(M.n_rows),
      orderMap(pc.orderMap),
      catMap(pc.catMap),
      colNames(pc.colNames),
      colTypes(pc.colTypes),
      masked(masked),
      mask(mask),
      subBoxFound(false) {};

  ColWorker(const ColWorker& cw, RcppParallel::Split)
    : M(cw.M),
      y(cw.y),
      p(cw.p), minSup(cw.minSup), N(cw.N),
      orderMap(cw.orderMap),
      catMap(cw.catMap),
      colNames(cw.colNames),
      colTypes(cw.colTypes),
      masked(cw.masked),
      mask(cw.mask),
      subBoxFound(false) {};


  void operator()(size_t begin, size_t end);
  void join(const ColWorker& cw);
  SubBox findNumCandidate(const arma::uword& colId);
  SubBox findCatCandidate(const arma::uword& colId);
};

std::string double2string(const double& d);
double mean(const std::vector<arma::uword>& idx, const arma::vec& y);


#endif
