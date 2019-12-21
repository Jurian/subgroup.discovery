// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
#include "prim.hpp"
#include "preprocess.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::uvec peel(const PrimContainer& pc) {

  std::string* colNames = pc.colNames;
  ColumnType* colTypes = pc.colTypes;
  std::map<arma::uword, arma::uvec> orderMap = pc.orderMap;
  std::map<arma::uword, std::map<arma::uword, Category>> catMap = pc.catMap;
  arma::mat M = pc.M;
  arma::vec y = pc.y;
  double alpha = pc.alpha;
  double minSup = pc.minSup;
  bool verbose = pc.verbose;

  bool* mask = new bool[M.n_rows]{};
  double boxSize = 1;
  const double N = M.n_rows;
  uword masked = 0;

  bool subBoxFound;

  do {
    SubBox bestSubBox, newSubBox;
    subBoxFound = false;
    for(uword i = 0; i < M.n_cols; i++) {

      try {

        switch(colTypes[i]) {
        case NUMERIC: {

          NumColWorker nc = {
            orderMap.find(i)->second, M.col(i), y, i, alpha, minSup, N, masked, mask
          };

          newSubBox = nc.findCandidate();
          break;
        }
        case FACTOR: {

          CatColWorker cc = {
            catMap.find(i)->second, M.col(i), y, i, alpha, minSup, N, masked, mask
          };

          newSubBox = cc.findCandidate();
          break;
        }
        default: stop("Invalid type for column %s", colNames[i]);
        }

        if(!subBoxFound || newSubBox.compare(bestSubBox) == 1) {
          bestSubBox = newSubBox;
          subBoxFound = true;
        }

      } catch(int e) {
        // No subBox was found for this column
      }
    }

    if(subBoxFound) {

      boxSize = bestSubBox.keep.size() / N;

      // Add the new subbox to the mask
      for(uword i = 0; i < bestSubBox.remove.size(); i++) {
        mask[bestSubBox.remove[i]] = true;
      }

      masked += bestSubBox.remove.size();

      if(verbose) {
        Rcout << "Box ratio: " << boxSize << endl;
        Rcout << "Rows to keep: " << bestSubBox.keep.size() << endl;
        Rcout << "Rows to remove: " << bestSubBox.remove.size() << endl;
        Rcout << "Rows masked: " << masked-bestSubBox.remove.size() << " + " << bestSubBox.remove.size() << " = ";
        Rcout << masked << endl;
        Rcout << "Peel off: " << colNames[bestSubBox.col];
        switch(bestSubBox.side) {
        case NUM_LEFT: Rcout << " < " << bestSubBox.value; break;
        case NUM_RIGHT: Rcout << " > " << bestSubBox.value; break;
        case CATEGORY: Rcout << " == " << bestSubBox.value; break;
        }
        Rcout << endl << "Quality of remainder: " << bestSubBox.quality << endl;
        Rcout << "----------------------------------" << endl;
      }
    }

  } while (subBoxFound);

  // All the rows which have not been masked belong to the final subbox
  vector<uword> outV;
  for(uword i = 0; i < M.n_rows; i++) {
    if(!mask[i]) outV.push_back(i+1);
  }

  return uvec(outV);
}

arma::uvec prim(const DataFrame& df, const vec& y, const double& alpha, const double& minSup, const bool& verbose) {

  PrimContainer pc = preprocess(df, y, alpha, minSup, verbose);

  return peel(pc);
}







