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
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <boost/dynamic_bitset.hpp>
#include "quantile.h"
#include "mapreduce.h"

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace boost;

void SubBox::applyBox(const NumericMatrix& M, dynamic_bitset<>& mask) const {

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

bool SubBox::isBetterThan(const SubBox& cmp) const {

  bool q = this->quality > cmp.quality;
  bool e = this->quality == cmp.quality;
  bool s = this->remove.size() < cmp.remove.size();

  return q || (e && s);
}

List SubBox::toList() const {
  return List::create(
    _["column"] = this->col,
    _["value"] = this->value,
    _["quality"] = this->quality,
    _["support"] = this->support,
    _["type"] = this->type
  );
}

SubBox SubBox::fromList(List list) {
  vector<int> rm;
  SubBox box = {
    rm,
    list["column"],
    list["value"],
    list["quality"],
    list["support"],
    list["type"]
  };
  return box;
}

SubBox ColWorker::findNumCandidate(const int& colId) {

  const RMatrix<double>::Column col = M.column(colId);
  const int N = M.nrow();

  const double leftQuantile = quantile7(col, RVector<int>(colOrders.find(colId)->second), alpha, N, masked, mask);
  const double rightQuantile = quantile7(col, RVector<int>(colOrders.find(colId)->second), 1-alpha, N, masked, mask);

  vector<int> leftRemove, rightRemove;
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

  SubBox left, right;
  bool leftFound = false, rightFound = false;
  const double leftSupport = leftKeepCount / (double)N;

  if(leftRemoveCount > 0 && leftSupport >= minSup) {
    leftFound = true;
    left = {
      leftRemove,
      colId,
      leftQuantile,
      leftMean,
      leftSupport,
      BOX_NUM_LEFT
    };
  }

  const double rightSupport = rightKeepCount / (double)N;
  if(rightRemoveCount > 0 && rightSupport >= minSup) {
    rightFound = true;
    right = {
      rightRemove,
      colId,
      rightQuantile,
      rightMean,
      rightSupport,
      BOX_NUM_RIGHT
    };
  }

  if(!leftFound && !rightFound) throw 0;
  if(!leftFound && rightFound) return right;
  if(leftFound && !rightFound) return left;
  return left.isBetterThan(right) ? left : right;
}

SubBox ColWorker::findCatCandidate(const int& colId) {

  const RMatrix<double>::Column col = M.column(colId);
  const int nCats = colCats.find(colId)->second;
  const int N = M.nrow();

  SubBox bestSubBox, newSubBox;
  bool boxFound = false;

  // Calculate the quality of each category
  for(int c = 0; c < nCats; c++) {

    double quality = 0;

    vector<int> remove;

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

      newSubBox = {
        remove,
        colId,
        (double) c,
        quality,
        support,
        BOX_CATEGORY
      };

      if(!boxFound || newSubBox.isBetterThan(bestSubBox)) {
        bestSubBox = newSubBox;
        boxFound = true;
      }
    }
  }

  if(!boxFound) throw 0;
  return bestSubBox;
}

void ColWorker::operator()(size_t begin, size_t end) {

  SubBox newSubBox;

  for(size_t i = begin; i < end; i++) {

    try {

      if(colTypes[i] == COL_NUMERIC) {
        newSubBox = findNumCandidate(i);
      }
      else if(colTypes[i] == COL_CATEGORICAL) {
        newSubBox = findCatCandidate(i);
      }
      else {
        stop("Invalid type for column %i", i+1);
      }

      if(!subBoxFound || newSubBox.isBetterThan(bestSubBox)) {

        bestSubBox = newSubBox;
        subBoxFound = true;
      }

    } catch(int e) {
      // No subBox was found for this column
    }
  }

}

void ColWorker::join(const ColWorker& cw) {

  if(cw.subBoxFound && (!subBoxFound || cw.bestSubBox.isBetterThan(bestSubBox))) {
    bestSubBox = cw.bestSubBox;
    subBoxFound = true;
  }
}
