#ifndef _QUANTILES_HPP
#define _QUANTILES_HPP

// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

using dCol = RMatrix<double>::Column;
using iVec = RVector<int>;

namespace quantile {
  double q(const double& gamma, const double& i, const double& j);

  double g(const int& n, const double& p, const double& m, const int& j);

  int j(const int& n, const double& p, const double& m);
}

int findMaskedIndex(const int& j, const int& N, const bool* mask);

double quantile7(const dCol& col, const iVec& order, const double& p, const int& N, const int& masked, const bool* mask);

double quantile2(const dCol& col, const iVec& order, const double& p, const int& N, const int& masked, const bool* mask);

#endif
