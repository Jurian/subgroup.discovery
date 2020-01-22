#ifndef _PRIM_HPP
#define _PRIM_HPP

// [[Rcpp::depends(RcppParallel)]]

#include <vector>
#include <map>
#include <Rcpp.h>
#include <RcppParallel.h>
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

IntegerVector sortIndex(const dCol& col);
int countCategories(const dCol& col);

vector<SubBox> findSubBoxes(
    const dMat& M,
    const dVec& y,
    const iVec& colTypes,
    const iMap& colCats,
    const vMap& colOrders,
    const double& alpha,
    const double& minSup);


//' PRIM Peel
//'
//' This function iteratively peels away a small portion of the matrix M
//' such that the mean of the remaining values in y is maximized.
//'
//' @param M The data to peel away from
//' @param y The labels to evaluate peels
//' @param colTypes Indicates which columns are numeric (0) and which are categorical (1)
//' @param alpha The peeling quantile
//' @param minSup The minimum allowed size of the remainder after a peel
//' @return A list of peeling steps
//' @author Jurian Baas
//  [[Rcpp::export]]
List peelCpp(
    const NumericMatrix& M,
    const NumericVector& y,
    const IntegerVector& colTypes,
    const double& alpha,
    const double& minSup);

//' PRIM Validate
//'
//' This function iteratively goes through the steps from the peeling process
//' and removes any trailing peels which do not improve the overall quality.
//'
//' @param peelSteps Peeling result from calling peel()
//' @param M The data to peel away from
//' @param y The labels to evaluate peels
//' @return A list of peeling steps
//' @author Jurian Baas
//  [[Rcpp::export]]
List predictCpp(
    const List& peelSteps,
    const NumericMatrix& M,
    const NumericVector& y);

//' PRIM Simplify Rules
//'
//' This function will go through all boxes that were found and tries to remove redundant ones,
//' as well as grouping them by column
//'
//' @param peelSteps A list of peeling steps
//' @param colTypes Indicates which columns are numeric (0) and which are categorical (1)
//' @param bestBoxIndex The (zero) index of the best box
//' @return A list of peeling steps
//' @author Jurian Baas
//  [[Rcpp::export]]
List simplifyRules(
    const List& peelSteps,
    const IntegerVector& colTypes,
    const int& bestBoxIndex);

#endif
