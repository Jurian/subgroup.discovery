#ifndef _COLWORKERS_HPP
#define _COLWORKERS_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <RcppArmadillo.h>
#include "quantile.hpp"

enum BoxType { NUM_LEFT, NUM_RIGHT, CATEGORY };

std::string double2string(const double& d);
double mean(const std::vector<arma::uword>& idx, const arma::vec& y);

struct SubBox {
  std::vector<arma::uword> remove;
  std::vector<arma::uword> keep;
  arma::uword col;
  std::string value;
  double quality;
  BoxType side;

  int compare(SubBox& cmp);
};

struct Category {
  arma::uword factor;
  std::string level;
  bool* index;
};

struct ColWorker {

  const arma::vec col, y;
  const arma::uword colId, masked;
  const double p, minSup, N;
  const bool* mask;

  ColWorker (
    const arma::vec& col, const arma::vec& y,
    const arma::uword& colId, const double& p, const double& minSup,
    const double& N, const arma::uword& masked, const bool* mask)
    : col(col), y(y), colId(colId), masked(masked), p(p), minSup(minSup), N(N), mask(mask) {};
};

struct CatColWorker : ColWorker {

  const std::map<arma::uword, Category> catMap;

  CatColWorker (
    const std::map<arma::uword, Category>& catMap,
    const arma::vec& col, const arma::vec& y,
    const arma::uword& colId, const double& p, const double& minSup,
    const double& N, const arma::uword& masked, const bool* mask
  ) : ColWorker (
    col, y, colId, p, minSup, N, masked, mask
  ), catMap(catMap) {}

  SubBox findCandidate();
};

struct NumColWorker : ColWorker {

  const arma::uvec order;

  NumColWorker(
    const arma::uvec& order,
    const arma::vec& col, const arma::vec& y,
    const arma::uword& colId, const double& p, const double& minSup,
    const double& N, const arma::uword& masked, const bool* mask
  ) : ColWorker (
      col, y, colId, p, minSup, N, masked, mask
  ), order(order) {}

  SubBox findCandidate();
};

#endif
