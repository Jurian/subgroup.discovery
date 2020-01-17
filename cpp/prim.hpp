#ifndef _PRIM_HPP
#define _PRIM_HPP

// [[Rcpp::depends(RcppArmadillo)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
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

IntegerVector sortIndex(const dCol& col);
int countCategories(const dCol& col);

vector<SubBox> findSubBoxes(
    const dMat& M,
    const dVec& y,
    const iVec& colTypes,
    const iMap& colCats,
    const vMap& colOrders,
    const double& alpha,
    const double& minSup,
    const bool& verbose);

// [[Rcpp::export]]
void peel(
    const NumericMatrix& M,
    const NumericVector& y,
    const IntegerVector& colTypes,
    const double& alpha,
    const double& minSup,
    const bool& verbose = false);


#endif
