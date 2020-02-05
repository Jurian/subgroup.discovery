/*
 * Subgroup Discovery
 * Copyright (C) 2020  Jurian Baas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/dynamic_bitset.hpp>
#include "quantile.h"

using namespace Rcpp;

double q(const double& gamma, const double& i, const double& j) {
  return (1 - gamma) * i + gamma * j;
}

double g(const int& n, const double& p, const double& m, const int& j) {
  return n * p + m - j;
}

int findMaskedIndex(const int& j, const int& N, const boost::dynamic_bitset<>& mask) {
  int k = 0;
  for(int i = 0; i < N; i++) {
    if(!mask[i]) k++;
    if(k == j) return i;
  }
  throw 0;
}

double quantile7 (const RcppParallel::RMatrix<double>::Column& col, const RcppParallel::RVector<int>& order, const double& p,
                  const int& N, const int& masked, const boost::dynamic_bitset<>& mask) {

  const int n = N - masked;

  const double index = 1 + (n - 1) * p;
  const int lo = floor(index);
  const int hi = ceil(index);
  const double gamma = index - lo;

  return q(gamma, col[order[findMaskedIndex(lo, N, mask)]], col[order[findMaskedIndex(hi, N, mask)]]);
}


double quantile2 (const RcppParallel::RMatrix<double>::Column& col, const RcppParallel::RVector<int>& order, const double& p,
                  const int& N, const int& masked, const boost::dynamic_bitset<>& mask) {

  // m is a constant determined by quantile type
  const double m = 0;
  const int n = N - masked;
  const int j = floor(n * p + m);
  const double gamma = g(n,p,m,j) == 0 ? 0.5 : 1;

  return q(gamma, col[order[findMaskedIndex(j, N, mask)]], col[order[findMaskedIndex(j + 1, N, mask)]]);
}

