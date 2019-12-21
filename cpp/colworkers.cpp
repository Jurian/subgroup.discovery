
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <RcppArmadillo.h>
#include "colworkers.hpp"
#include "quantile.hpp"

static bool verbose = false;

std::string double2string(const double& d) {
  std::ostringstream strs;
  strs << d;
  return strs.str();
}

double mean(const std::vector<arma::uword>& idx, const arma::vec& y) {
  double sum = 0;
  for(std::vector<arma::uword>::const_iterator it = idx.begin(); it != idx.end(); ++it) {
    sum += y[*it];
  }
  return sum/idx.size();
}

int SubBox::compare(SubBox& cmp) {

  bool q = quality > cmp.quality;
  bool e = quality == cmp.quality;
  bool s = keep.size() > cmp.keep.size();

  if(q || (e && s)) return 1;
  else return -1;
}

SubBox CatColWorker::findCandidate() {

  SubBox bestSubBox, newSubBox;
  bool subBoxFound = false;

  // Calculate the quality of each category
  std::map<arma::uword, Category>::const_iterator it;
  for(it = catMap.begin(); it != catMap.end(); ++it) {

    Category c = it->second;
    double quality = 0;

    std::vector<arma::uword> keep, remove;

    arma::uword k = 0;
    for(arma::uword i = 0; i < col.n_elem; i++) {
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

  if(verbose) {
    if(subBoxFound) {
      Rcpp::Rcout << "Found categorical subbox:" << std::endl;
      Rcpp::Rcout << "Column #" << bestSubBox.col << ", remove category " << bestSubBox.value << std::endl;
      Rcpp::Rcout << "New quality: " << bestSubBox.quality << std::endl << std::endl;
    } else {
      Rcpp::Rcout << "No valid categorical subbox for column #" << colId << std::endl;
    }
  }

  if(!subBoxFound) throw 0;
  return bestSubBox;
};

SubBox NumColWorker::findCandidate() {

  const double leftValue = quantile7(col, order, p, N, masked, mask);
  const double rightValue = quantile7(col, order, 1-p, N, masked, mask);

  std::vector<arma::uword> leftRemove, leftKeep, rightRemove, rightKeep;

  for(arma::uword i = 0; i < N; i++) {
    if(mask[i]) continue;
    if(col[i] < leftValue) leftRemove.push_back(i);
    else leftKeep.push_back(i);
    if(col[i] > rightValue) rightRemove.push_back(i);
    else rightKeep.push_back(i);
  }

  const arma::uword leftRemoveCount = leftRemove.size();
  const arma::uword leftKeepCount = leftKeep.size();
  const arma::uword rightRemoveCount = rightRemove.size();
  const arma::uword rightKeepCount = rightKeep.size();

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

  if(verbose) {
    if(leftFound) {
      Rcpp::Rcout << "Found left subbox:" << std::endl;
      Rcpp::Rcout << "Column #" << colId << ", remove values < " << left.value << std::endl;
      Rcpp::Rcout << "New quality: " << left.quality << std::endl << std::endl;
    }

    if(rightFound) {
      Rcpp::Rcout << "Found right subbox:" << std::endl;
      Rcpp::Rcout << "Column #" << colId << ", remove values > " << right.value << std::endl;
      Rcpp::Rcout << "New quality: " << right.quality << std::endl << std::endl;
    }
  }

  if(!leftFound && !rightFound) throw 0;
  return left.compare(right) == 1 ? left : right;
};
