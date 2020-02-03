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

#ifndef _QUANTILES_HPP
#define _QUANTILES_HPP

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/dynamic_bitset.hpp>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace boost;

double q(const double& gamma, const double& i, const double& j);

double g(const int& n, const double& p, const double& m, const int& j);

int findMaskedIndex(const int& j, const int& N, const dynamic_bitset<>& mask);

double quantile7(const RMatrix<double>::Column& col, const RVector<int>& order, const double& p, const int& N, const int& masked, const dynamic_bitset<>& mask);

double quantile2(const RMatrix<double>::Column& col, const RVector<int>& order, const double& p, const int& N, const int& masked, const dynamic_bitset<>& mask);

#endif
