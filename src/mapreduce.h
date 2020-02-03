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

#ifndef _MAPREDUCE_HPP
#define _MAPREDUCE_HPP

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/dynamic_bitset.hpp>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace boost;

static const int COL_NUMERIC = 0;
static const int COL_CATEGORICAL = 1;

static const int BOX_NUM_LEFT = 0;
static const int BOX_NUM_RIGHT = 1;
static const int BOX_CATEGORY = 2;

struct SubBox {
  vector<int> remove;
  int col;
  double value;
  double quality;
  double support;
  int type;

  void applyBox(const NumericMatrix& M, dynamic_bitset<>& mask) const;
  bool isBetterThan(const SubBox& cmp) const;
  List toList() const;
  static SubBox fromList(List list);
};

struct ColWorker : public Worker {

  const RMatrix<double>& M;
  const RVector<double>& y;
  const RVector<int>& colTypes;
  const map<int, int>& colCats;
  const map<int, IntegerVector>& colOrders;
  const double& alpha;
  const double& minSup;
  const int& masked;
  const dynamic_bitset<>& mask;

  bool subBoxFound;
  SubBox bestSubBox;

  // constructors
  ColWorker (
      const RMatrix<double>& M,
      const RVector<double>& y,
      const RVector<int>& colTypes,
      const map<int, int>& colCats,
      const map<int, IntegerVector>& colOrders,
      const double& alpha,
      const double& minSup,
      const int& masked,
      const dynamic_bitset<>& mask)
    : M(M),
      y(y),
      colTypes(colTypes),
      colCats(colCats),
      colOrders(colOrders),
      alpha(alpha),
      minSup(minSup),
      masked(masked),
      mask(mask),
      subBoxFound(false) {};

  ColWorker(const ColWorker& cw, RcppParallel::Split)
    : M(cw.M),
      y(cw.y),
      colTypes(cw.colTypes),
      colCats(cw.colCats),
      colOrders(cw.colOrders),
      alpha(cw.alpha),
      minSup(cw.minSup),
      masked(cw.masked),
      mask(cw.mask),
      subBoxFound(cw.subBoxFound),
      bestSubBox(cw.bestSubBox) {};

  void operator()(size_t begin, size_t end);
  void join(const ColWorker& cw);
  SubBox findNumCandidate(const int& colId);
  SubBox findCatCandidate(const int& colId);
};

#endif
