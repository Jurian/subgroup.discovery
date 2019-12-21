#include "quantile.hpp"
#include <RcppArmadillo.h>

namespace quantile {
  double q(const double& gamma, const double& i, const double& j) {
    return (1 - gamma) * i + gamma * j;
  }

  double g(const arma::uword& n, const double& p, const double& m, const arma::uword& j) {
    return n * p + m - j;
  }

  arma::uword j(const arma::uword& n, const double& p, const double& m) {
    return floor(n * p + m);
  }
}

arma::uword findMaskedIndex(const arma::uword& j, const arma::uword& N, const bool* mask) {
  arma::uword k = 0;
  for(arma::uword i = 0; i < N; i++) {
    if(!mask[i]) k++;
    if(k == j) return i;
  }
  throw 0;
}

double quantile7 (const arma::vec& col, const arma::uvec& order, const double& p,
                  const arma::uword& N,const arma::uword& masked, const bool* mask) {

  const arma::uword n = N - masked;

  const double index = 1 + (n - 1) * p;
  const arma::uword lo = floor(index);
  const arma::uword hi = ceil(index);
  const double gamma = index - lo;

  return quantile::q(gamma, col[order[findMaskedIndex(lo, N, mask)]], col[order[findMaskedIndex(hi, N, mask)]]);
}


double quantile2 (const arma::vec& col, const arma::uvec& order, const double& p,
                  const arma::uword& N, const arma::uword& masked, const bool* mask) {

  // m is a constant determined by quantile type
  const double m = 0;
  const arma::uword n = N - masked;
  const arma::uword j = quantile::j(n,p,m);
  const double gamma = quantile::g(n,p,m,j) == 0 ? 0.5 : 1;

  return quantile::q(gamma, col[order[findMaskedIndex(j, N, mask)]], col[order[findMaskedIndex(j + 1, N, mask)]]);
}

