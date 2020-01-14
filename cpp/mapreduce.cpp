// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "quantile.hpp"
#include "mapreduce.hpp"

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;
using namespace arma;

string double2string(const double& d) {
  ostringstream strs;
  strs << d;
  return strs.str();
}

double mean(const vector<uword>& idx, const vec& y) {
  double sum = 0;
  for(vector<uword>::const_iterator it = idx.begin(); it != idx.end(); ++it) {
    sum += y[*it];
  }
  return sum/idx.size();
}

int SubBox::compare(const SubBox& cmp) const {

  bool q = quality > cmp.quality;
  bool e = quality == cmp.quality;
  bool s = keep.size() > cmp.keep.size();

  if(q || (e && s)) return 1;
  else return -1;
}

SubBox ColWorker::findNumCandidate(const uword& colId) {

  const vec col = M.col(colId);

  const double leftValue = quantile7(col, orderMap.find(colId)->second, p, N, masked, mask);
  const double rightValue = quantile7(col, orderMap.find(colId)->second, 1-p, N, masked, mask);

  vector<uword> leftRemove, leftKeep, rightRemove, rightKeep;

  for(uword i = 0; i < N; i++) {
    if(mask[i]) continue;
    if(col[i] < leftValue) leftRemove.push_back(i);
    else leftKeep.push_back(i);
    if(col[i] > rightValue) rightRemove.push_back(i);
    else rightKeep.push_back(i);
  }

  const uword leftRemoveCount = leftRemove.size();
  const uword leftKeepCount = leftKeep.size();
  const uword rightRemoveCount = rightRemove.size();
  const uword rightKeepCount = rightKeep.size();

  SubBox left, right;
  bool leftFound = false, rightFound = false;

  if(leftRemoveCount > 0 && leftKeepCount / N >= minSup) {
    leftFound = true;
    left = {
      leftRemove,
      leftKeep,
      colId,
      double2string(leftValue),
      mean(leftKeep, y),
      NUM_LEFT
    };
  }

  if(rightRemoveCount > 0 && rightKeepCount / N >= minSup) {
    rightFound = true;
    right = {
      rightRemove,
      rightKeep,
      colId,
      double2string(rightValue),
      mean(rightKeep, y),
      NUM_RIGHT
    };
  }

  if(!leftFound && !rightFound) throw 0;
  return left.compare(right) == 1 ? left : right;
};

SubBox ColWorker::findCatCandidate(const uword& colId) {

  const vec col = M.col(colId);
  const map<uword, Category> categories = catMap.find(colId)->second;

  SubBox bestSubBox, newSubBox;
  bool subBoxFound = false;

  // Calculate the quality of each category
  map<uword, Category>::const_iterator it;
  for(it = categories.begin(); it != categories.end(); ++it) {

    Category c = it->second;
    double quality = 0;

    vector<uword> keep, remove;

    uword k = 0;
    for(uword i = 0; i < col.n_elem; i++) {
      if(mask[i]) continue; // Skip masked rows
      if(!c.index[i]) { // Take cumulative moving avarage of all other categories
        quality = quality + ((y[i] - quality) / (k + 1));
        keep.push_back(i);
        k++;
      }else {
        remove.push_back(i);
      }
    }

    if(keep.size() == 0) continue;
    if(remove.size() == 0) continue;

    const double removeFraction = remove.size() / (double) (remove.size() + keep.size());
    if(removeFraction < p && keep.size() / N >= minSup) {

      newSubBox = {
        remove,
        keep,
        colId,
        c.level,
        quality,
        CATEGORY
      };

      if(!subBoxFound || newSubBox.compare(bestSubBox) == 1) {
        bestSubBox = newSubBox;
        subBoxFound = true;
      };
    }
  }

  if(!subBoxFound) throw 0;
  return bestSubBox;
};

void ColWorker::operator()(size_t begin, size_t end) {

  SubBox newSubBox;

  for(size_t i = begin; i < end; i++) {

    try {

      if(colTypes[i] == NUMERIC) {
        newSubBox = findNumCandidate(i);
      }
      else if(colTypes[i] == FACTOR) {
        newSubBox = findCatCandidate(i);
      }
      else {
        stop("Invalid type for column %s", colNames[i]);
      }

      if(!subBoxFound || newSubBox.compare(bestSubBox) == 1) {
        bestSubBox = newSubBox;
        subBoxFound = true;
      }

    } catch(int e) {
      // No subBox was found for this column
    }
  }
};

void ColWorker::join(const ColWorker& cw) {

  if(cw.subBoxFound && (!subBoxFound || cw.bestSubBox.compare(bestSubBox) == 1)) {
    bestSubBox = cw.bestSubBox;
    subBoxFound = true;
  }
};
