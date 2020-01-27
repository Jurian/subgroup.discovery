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

IntegerVector sortIndex(const NumericMatrix::ConstColumn& col);
int countCategories(const NumericMatrix::ConstColumn& col);

List findSubBoxes(
    const RMatrix<double>& M,
    const RVector<double>& y,
    const RVector<int>& colTypes,
    const map<int, int>& colCats,
    const map<int, IntegerVector>& colOrders,
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
//' This function evaluates all the steps from the peeling process on new data.
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
