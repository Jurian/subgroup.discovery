// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "prim.hpp"
#include "mapreduce.hpp"

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
  //return iVec(index);
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
    const double& minSup,
    const bool& verbose) {

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

      if(verbose) {
        int keepSize = M.nrow() - masked;
        double boxSize = keepSize / M.nrow() ;

        Rcout << M.nrow() << endl;
        Rcout << masked << endl;
        Rcout << keepSize << endl;
        Rcout << boxSize << endl << endl;
      }


      peelSteps.push_back(bestSubBox);
    }

  } while (subBoxFound);

  return peelSteps;
}

void peel (
    const NumericMatrix& M,
    const NumericVector& y,
    const IntegerVector& colTypes,
    const double& alpha,
    const double& minSup,
    const bool& verbose) {

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

  vector<SubBox> peelSteps = findSubBoxes(
    dMat(M),
    dVec(y),
    iVec(colTypes),
    colCats,
    colOrders,
    alpha,
    minSup,
    verbose
  );
  //return peel(pc);
}

