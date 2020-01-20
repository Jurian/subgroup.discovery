#ifndef _MAPREDUCE_HPP
#define _MAPREDUCE_HPP

// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

using dMat = RMatrix<double>;
using dVec = RVector<double>;
using iVec = RVector<int>;
using iMap = map<int, int>;
using vMap = map<int, IntegerVector>;

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

  bool isBetterThan(const SubBox& cmp) const;
  List toList() const;
};

struct ColWorker : public Worker {

  const dMat& M;
  const dVec& y;
  const iVec& colTypes;
  const iMap& colCats;
  const vMap& colOrders;
  const double& alpha;
  const double& minSup;
  const int& masked;
  const bool* mask;

  bool subBoxFound;
  SubBox bestSubBox;

  // constructors
  ColWorker (
      const dMat& M,
      const dVec& y,
      const iVec& colTypes,
      const iMap& colCats,
      const vMap& colOrders,
      const double& alpha,
      const double& minSup,
      const int& masked,
      const bool* mask)
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
