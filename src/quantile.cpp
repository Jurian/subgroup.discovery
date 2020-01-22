// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include "quantile.h"

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

using dCol = RMatrix<double>::Column;
using iVec = RVector<int>;

namespace quantile {
  double q(const double& gamma, const double& i, const double& j) {
    return (1 - gamma) * i + gamma * j;
  }

  double g(const int& n, const double& p, const double& m, const int& j) {
    return n * p + m - j;
  }

  int j(const int& n, const double& p, const double& m) {
    return floor(n * p + m);
  }
}

int findMaskedIndex(const int& j, const int& N, const bool* mask) {
  int k = 0;
  for(int i = 0; i < N; i++) {
    if(!mask[i]) k++;
    if(k == j) return i;
  }
  throw 0;
}

double quantile7 (const dCol& col, const iVec& order, const double& p,
                  const int& N,const int& masked, const bool* mask) {

  const int n = N - masked;

  const double index = 1 + (n - 1) * p;
  const int lo = floor(index);
  const int hi = ceil(index);
  const double gamma = index - lo;

  return quantile::q(gamma, col[order[findMaskedIndex(lo, N, mask)]], col[order[findMaskedIndex(hi, N, mask)]]);
}


double quantile2 (const dCol& col, const iVec& order, const double& p,
                  const int& N, const int& masked, const bool* mask) {

  // m is a constant determined by quantile type
  const double m = 0;
  const int n = N - masked;
  const int j = quantile::j(n,p,m);
  const double gamma = quantile::g(n,p,m,j) == 0 ? 0.5 : 1;

  return quantile::q(gamma, col[order[findMaskedIndex(j, N, mask)]], col[order[findMaskedIndex(j + 1, N, mask)]]);
}

