// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "prim.h"
#include "mapreduce.h"

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

using dMat = RMatrix<double>;
using dVec = RVector<double>;
using iVec = RVector<int>;
using iMap = map<int, int>;
using vMap = map<int, IntegerVector>;
using dCol = NumericMatrix::ConstColumn;

IntegerVector sortIndex(const dCol& col) {
  IntegerVector index(col.size());
  for (int i = 0 ; i != index.size() ; i++) {
    index[i] = i;
  }
  sort(index.begin(), index.end(),
       [&](const int& a, const int& b) {
         return (col[a] < col[b]);
       }
  );
  return index;
}

int countCategories(const dCol& col) {
  return unique(col).length();
}

vector<SubBox> findSubBoxes(
    const dMat& M,
    const dVec& y,
    const iVec& colTypes,
    const iMap& colCats,
    const vMap& colOrders,
    const double& alpha,
    const double& minSup) {

  bool* mask = new bool[M.nrow()]{};
  int masked = 0;

  bool subBoxFound;

  vector<SubBox> peelSteps;

  do {

    ColWorker cw(
        M,
        y,
        colTypes,
        colCats,
        colOrders,
        alpha,
        minSup,
        masked,
        mask);

    parallelReduce(0, M.ncol(), cw);

    subBoxFound = cw.subBoxFound;

    if(subBoxFound) {
      SubBox bestSubBox = cw.bestSubBox;

      // Add the new subbox to the mask
      for(size_t i = 0; i < bestSubBox.remove.size(); i++) {
        mask[bestSubBox.remove[i]] = true;
      }

      masked += bestSubBox.remove.size();
      peelSteps.push_back(bestSubBox);
    }

  } while (subBoxFound);

  return peelSteps;
}

List peelCpp (
    const NumericMatrix& M,
    const NumericVector& y,
    const IntegerVector& colTypes,
    const double& alpha,
    const double& minSup) {

  iMap colCats;
  vMap colOrders;

  for(int i = 0; i < colTypes.length(); i++ ) {

    if(colTypes[i] == COL_NUMERIC) {

      colOrders[i] = sortIndex(M(_, i));

    } else if(colTypes[i] == COL_CATEGORICAL) {

      colCats[i] = countCategories(M(_, i));

    } else {
      stop("Column nr %i is not numeric or categorical", i);
    }
  }

  vector<SubBox> boxes = findSubBoxes(
    dMat(M),
    dVec(y),
    iVec(colTypes),
    colCats,
    colOrders,
    alpha,
    minSup
  );

  List peelSteps = List::create();
  for(vector<SubBox>::iterator it = boxes.begin(); it != boxes.end(); ++it) {
    peelSteps.push_back((*it).toList());
  }

  return peelSteps;
}

List predictCpp (
    const List& peelSteps,
    const NumericMatrix& M,
    const NumericVector& y) {

  const int N = M.nrow();

  int masked = 0;
  bool* mask = new bool[M.nrow()]{};

  List validationSteps = List::create();

  for(List::const_iterator it = peelSteps.begin(); it != peelSteps.end(); ++it) {

    const List peel = *it;
    const int colId = peel["column"];
    const int type = peel["type"];
    const double value = peel["value"];

    double quality = 0;

    int k = 1;
    vector<int> remove;

    dCol column = M( _, colId);

    if(type == BOX_NUM_LEFT) {

      for(int i = 0; i < column.size(); i++) {

        if(mask[i]) continue;

        if(column[i] < value) {
          remove.push_back(i);
        } else {
          // Cumulative moving avarage of rows we keep
          quality += (y[i] - quality) / k++;
        }

      }

    } else if (type == BOX_NUM_RIGHT) {

      for(int i = 0; i < column.size(); i++) {

        if(mask[i]) continue;

        if(column[i] > value) {
          remove.push_back(i);
        } else {
          // Cumulative moving avarage of rows we keep
          quality += (y[i] - quality) / k++;
        }

      }

    } else if (type == BOX_CATEGORY) {

      for(int i = 0; i < column.size(); i++) {

        if(mask[i]) continue;

        if(column[i] == value) {
          remove.push_back(i);
        } else {
          // Cumulative moving avarage of rows we keep
          quality += (y[i] - quality) / k++;
        }

      }

    }

    // Add the new subbox to the mask
    for(size_t i = 0; i < remove.size(); i++) {
      mask[remove[i]] = true;
    }
    masked += remove.size();

    const int keepCount = N - masked;
    const double support = keepCount / (double)N;

    validationSteps.push_back(
      List::create(
        _["column"] = colId,
        _["type"] = type,
        _["value"] = value,
        _["quality"] = quality,
        _["support"] = support
      )
    );

  }

  IntegerVector finalBoxIndex;
  for(int i = 0; i < N; i++) {
    if(mask[i]) continue;
    finalBoxIndex.push_back(i);
  }

  List validationResult = List::create();

  validationResult.push_back(validationSteps);
  validationResult.push_back(finalBoxIndex);

  return validationResult;

}

