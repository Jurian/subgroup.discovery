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

#include <vector>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/dynamic_bitset.hpp>
#include "quantile.h"
#include "mapreduce.h"

using namespace Rcpp;

void Peel::apply(const NumericMatrix& M, boost::dynamic_bitset<>& mask) const {

  const size_t N = M.nrow();
  NumericMatrix::ConstColumn column = M(_, this->col);

  if(this->type == BOX_NUM_LEFT) {
    for(size_t i = 0; i < N; i++) {
      if(column[i] < this->value) mask.set(i);
    }
  }

  else if(this->type == BOX_NUM_RIGHT) {
    for(size_t i = 0; i < N; i++) {
      if(column[i] > this->value) mask.set(i);
    }
  }

  else if(this->type == BOX_CATEGORY) {
    for(size_t i = 0; i < N; i++) {
      if(column[i] == this->value) mask.set(i);
    }
  }
}

bool Peel::isBetterThan(const Peel& cmp) const {

  // When we choose between two peels we always choose
  // a valid one over an invalid one. If both are invalid,
  // then the peel is not "better".
  // If both are valid, choose the peel with the highest
  // quality and break ties on support

  if(this->valid && !cmp.valid) return true;
  if(!this->valid && cmp.valid) return false;
  if(!this->valid && !cmp.valid) return false;

  bool q = this->quality > cmp.quality;
  bool e = this->quality == cmp.quality;
  bool s = this->remove.size() < cmp.remove.size();

  return q || (e && s);
}

List Peel::toList() const {
  return List::create(
    _["column"] = this->col,
    _["value"] = this->value,
    _["quality"] = this->quality,
    _["support"] = this->support,
    _["type"] = this->type,
    _["valid"] = this->valid
  );
}

Peel Peel::fromList(List list) {

  Peel box = {
    (bool)list["valid"],
    (int)list["column"],
    (int)list["type"],
    (double)list["value"],
    (double)list["quality"],
    (double)list["support"],
    std::vector<int>()
  };
  return box;
}

Peel ColWorker::findNumCandidate(const int& colId) {

  const RcppParallel::RMatrix<double>::Column col = M.column(colId);
  const int N = M.nrow();

  const double leftQuantile = quantile7(col, RcppParallel::RVector<int>(colOrders.find(colId)->second), alpha, N, masked, mask);
  const double rightQuantile = quantile7(col, RcppParallel::RVector<int>(colOrders.find(colId)->second), 1-alpha, N, masked, mask);

  std::vector<int> leftRemove, rightRemove;
  double leftMean = 0, rightMean = 0;
  int leftK = 1, rightK = 1;

  for(int i = 0; i < N; i++) {
    if(mask[i]) continue;

    if(col[i] < leftQuantile) {
      leftRemove.push_back(i);
    } else {
      // Cumulative moving avarage of rows we keep
      leftMean += (y[i] - leftMean) / leftK++;
    }

    if(col[i] > rightQuantile) {
      rightRemove.push_back(i);
    } else {
      // Cumulative moving avarage of rows we keep
      rightMean += (y[i] - rightMean) / rightK++;
    }
  }

  const int leftRemoveCount = leftRemove.size();
  const int leftKeepCount = N - masked - leftRemove.size();
  const int rightRemoveCount = rightRemove.size();
  const int rightKeepCount = N - masked - rightRemove.size();

  // Create empty invalid boxes for both left and right
  Peel left = {
    colId,
    BOX_NUM_LEFT
  };

  Peel right = {
    colId,
    BOX_NUM_RIGHT
  };

  const double leftSupport = leftKeepCount / (double)N;
  if(leftRemoveCount > 0 && leftSupport >= minSup) {
    left.valid = true;
    left.value = leftQuantile;
    left.quality = leftMean;
    left.support = leftSupport;
    left.remove = leftRemove;
  }

  const double rightSupport = rightKeepCount / (double)N;
  if(rightRemoveCount > 0 && rightSupport >= minSup) {
    right.valid = true;
    right.value = rightQuantile;
    right.quality = rightMean;
    right.support = rightSupport;
    right.remove = rightRemove;
  }

  // Note that we might still return an invalid box here
  return left.isBetterThan(right) ? left : right;
}

Peel ColWorker::findCatCandidate(const int& colId) {

  const RcppParallel::RMatrix<double>::Column col = M.column(colId);
  const int nCats = colCats.find(colId)->second;
  const int N = M.nrow();

  // Create a default invalid box
  Peel box = {
    colId,
    BOX_CATEGORY
  };

  double bestSupport = 0;
  double bestQuality = 0;

  // Calculate the quality of each category
  for(int c = 0; c < nCats; c++) {

    double quality = 0;

    std::vector<int> remove;

    int k = 1;
    for(size_t i = 0; i < col.size(); i++) {
      if(mask[i]) continue; // Skip masked rows
      if((int) col[i] != c) {
        // Take cumulative moving avarage of all other categories
        quality += (y[i] - quality) / k++;
      }else {
        remove.push_back(i);
      }
    }

    const int keepSize = N - masked - remove.size();

    if(keepSize == 0) continue;
    if(remove.size() == 0) continue;

    const double removeFraction = remove.size() / (double) (remove.size() + keepSize);
    const double support = keepSize / (double)N;

    if(removeFraction < alpha && support >= minSup) {

      if(quality > bestQuality || (quality == bestQuality && support > bestSupport) ) {

        box.valid = true;
        box.value = (double) c;
        box.quality = quality;
        box.support = support;
        box.remove = remove;

        bestSupport = support;
        bestQuality = quality;
      }
    }
  }

  // Note that we might still return an invalid box here
  return box;
}

void ColWorker::operator()(size_t begin, size_t end) {

  // We start with an invalid default box
  for(size_t i = begin; i < end; i++) {

    if(colTypes[i] == COL_NUMERIC) {

      Peel newPeel = findNumCandidate(i);

      // At this point the new box might be invalid!
      // Note that a valid box is always "better" than an invalid box
      if(newPeel.valid && newPeel.isBetterThan(this->bestPeel)) {
        this->bestPeel = newPeel;
      }

    } else {

      Peel newPeel = findCatCandidate(i);

      // At this point the new box might be invalid!
      // Note that a valid box is always "better" than an invalid box
      if(newPeel.valid && newPeel.isBetterThan(this->bestPeel)) {
        this->bestPeel = newPeel;
      }
    }
  }
}

void ColWorker::join(const ColWorker& cw) {
  // Note that a valid box is always "better" than an invalid box
  if(cw.bestPeel.isBetterThan(this->bestPeel)) {
    this->bestPeel = cw.bestPeel;
  }
}
