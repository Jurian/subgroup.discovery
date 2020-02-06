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

using namespace Rcpp;

IntegerVector sortIndex(const NumericMatrix::ConstColumn& col);
int countCategories(const NumericMatrix::ConstColumn& col);

List findPeels(
    const RcppParallel::RMatrix<double>& M,
    const RcppParallel::RVector<double>& y,
    const RcppParallel::RVector<int>& colTypes,
    const std::map<int, int>& colCats,
    const std::map<int, IntegerVector>& colOrders,
    const double& alpha,
    const double& minSup);

//  [[Rcpp::export]]
List peelCpp(
    const NumericMatrix& M,
    const NumericVector& y,
    const IntegerVector& colTypes,
    const double& alpha,
    const double& minSup);

//  [[Rcpp::export]]
List predictCpp(
    const List& peelResult,
    const NumericMatrix& M,
    const NumericVector& y);

//  [[Rcpp::export]]
List simplifyCpp(
    const List& peelSteps,
    const IntegerVector& colTypes,
    const size_t& boxId);

//  [[Rcpp::export]]
IntegerVector indexCpp(
        const List& boxes,
        const NumericMatrix& M,
        const size_t& boxId);

#endif
