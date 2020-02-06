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

using namespace Rcpp;

static const int COL_NUMERIC = 0;
static const int COL_CATEGORICAL = 1;

static const int BOX_NUM_LEFT = 0;
static const int BOX_NUM_RIGHT = 1;
static const int BOX_CATEGORY = 2;

struct Peel {

  bool valid;
  int col;
  int type;
  double value;
  double quality;
  double support;
  std::vector<int> remove;

  Peel()
  : valid(false),
    col(0),
    type(0),
    value(0),
    quality(0),
    support(0),
    remove(std::vector<int>()) {};

  Peel(
    int col,
    int type)
  : valid(false),
    col(col),
    type(type),
    value(0),
    quality(0),
    support(0),
    remove(std::vector<int>()) {};

  Peel(
    bool valid,
    int col,
    int type,
    double value,
    double quality,
    double support,
    std::vector<int> remove)
  : valid(valid),
    col(col),
    type(type),
    value(value),
    quality(quality),
    support(support),
    remove(remove) {};

  void apply(const NumericMatrix& M, boost::dynamic_bitset<>& mask) const;
  bool isBetterThan(const Peel& cmp) const;
  List toList() const;
  static Peel fromList(List list);
};

struct ColWorker : public RcppParallel::Worker {

  const RcppParallel::RMatrix<double>& M;
  const RcppParallel::RVector<double>& y;
  const RcppParallel::RVector<int>& colTypes;
  const std::map<int, int>& colCats;
  const std::map<int, IntegerVector>& colOrders;
  const double& alpha;
  const double& minSup;
  const int& masked;
  const boost::dynamic_bitset<>& mask;

  Peel bestPeel = {};

  // constructors
  ColWorker (
      const RcppParallel::RMatrix<double>& M,
      const RcppParallel::RVector<double>& y,
      const RcppParallel::RVector<int>& colTypes,
      const std::map<int, int>& colCats,
      const std::map<int, IntegerVector>& colOrders,
      const double& alpha,
      const double& minSup,
      const int& masked,
      const boost::dynamic_bitset<>& mask)
    : M(M),
      y(y),
      colTypes(colTypes),
      colCats(colCats),
      colOrders(colOrders),
      alpha(alpha),
      minSup(minSup),
      masked(masked),
      mask(mask) {};

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
      bestPeel(cw.bestPeel) {};

  void operator()(size_t begin, size_t end);
  void join(const ColWorker& cw);
  Peel findNumCandidate(const int& colId);
  Peel findCatCandidate(const int& colId);
};

#endif
