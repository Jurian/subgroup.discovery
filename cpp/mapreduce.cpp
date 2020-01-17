// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "quantile.hpp"
#include "mapreduce.hpp"

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

using pCol = RMatrix<double>::Column;
using iVec = RVector<int>;

bool SubBox::isBetterThan(const SubBox& cmp) const {

  bool q = quality > cmp.quality;
  bool e = quality == cmp.quality;
  bool s = remove.size() < cmp.remove.size();

  return q || (e && s);
}

SubBox ColWorker::findNumCandidate(const int& colId) {

  const pCol col = M.column(colId);
  const int N = M.nrow();

  const double leftQuantile = quantile7(col, iVec(colOrders.find(colId)->second), alpha, N, masked, mask);
  const double rightQuantile = quantile7(col, iVec(colOrders.find(colId)->second), 1-alpha, N, masked, mask);

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

  if(leftRemoveCount > 0 && leftKeepCount / N >= minSup) {
    leftFound = true;
    left = {
      leftRemove,
      colId,
      leftQuantile,
      leftMean,
      BOX_NUM_LEFT
    };
  }

  if(rightRemoveCount > 0 && rightKeepCount / N >= minSup) {
    rightFound = true;
    right = {
      rightRemove,
      colId,
      rightQuantile,
      rightMean,
      BOX_NUM_RIGHT
    };
  }

  if(!leftFound && !rightFound) throw 0;
  if(!leftFound && rightFound) return right;
  if(leftFound && !rightFound) return left;
  return left.isBetterThan(right) ? left : right;
};

SubBox ColWorker::findCatCandidate(const int& colId) {

  const pCol col = M.column(colId);
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
    if(removeFraction < alpha && keepSize / N >= minSup) {

      newSubBox = {
        remove,
        colId,
        (double) c,
        quality,
        BOX_CATEGORY
      };

      if(!boxFound || newSubBox.isBetterThan(bestSubBox)) {
        bestSubBox = newSubBox;
        boxFound = true;
      };
    }
  }

  if(!boxFound) throw 0;
  return bestSubBox;
};

void ColWorker::operator()(size_t begin, size_t end) {

  SubBox newSubBox;
  subBoxFound = false;
  for(size_t i = begin; i < end; i++) {

    try {
      //Rcout << i;
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
};

void ColWorker::join(const ColWorker& cw) {

  if(cw.subBoxFound && (!subBoxFound || cw.bestSubBox.isBetterThan(bestSubBox))) {
    bestSubBox = cw.bestSubBox;
    subBoxFound = true;
  }
};
