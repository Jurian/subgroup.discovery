#ifndef _QUANTILES_HPP
#define _QUANTILES_HPP

#include <RcppArmadillo.h>

namespace quantile {
  double q(const double& gamma, const double& i, const double& j);

  double g(const arma::uword& n, const double& p, const double& m, const arma::uword& j);

  arma::uword j(const arma::uword& n, const double& p, const double& m);
}

arma::uword findMaskedIndex(const arma::uword& j, const arma::uword& N, const bool* mask);

double quantile7(const arma::vec& col, const arma::uvec& order, const double& p, const arma::uword& N, const arma::uword& masked, const bool* mask);

double quantile2(const arma::vec& col, const arma::uvec& order, const double& p, const arma::uword& N, const arma::uword& masked, const bool* mask);

#endif
