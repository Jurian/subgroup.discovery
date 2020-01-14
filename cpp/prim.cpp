// [[Rcpp::depends(RcppArmadillo)]]

#include <string>
#include <vector>
#include <map>
#include <RcppArmadillo.h>
#include "prim.hpp"
#include "preprocess.hpp"
#include "mapreduce.hpp"

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

arma::uvec peel(const PrimContainer& pc) {

  std::string* colNames = pc.colNames;
  //ColumnType* colTypes = pc.colTypes;
  //std::map<arma::uword, arma::uvec> orderMap = pc.orderMap;
  //std::map<arma::uword, std::map<arma::uword, Category>> catMap = pc.catMap;
  arma::mat M = pc.M;
  //arma::vec y = pc.y;
  //double alpha = pc.alpha;
  //double minSup = pc.minSup;
  bool verbose = pc.verbose;

  bool* mask = new bool[M.n_rows]{};
  double boxSize = 1;
  const double N = M.n_rows;
  uword masked = 0;

  bool subBoxFound;

  do {

    ColWorker cw(pc, masked, mask);

    parallelReduce(0, M.n_cols, cw);

    subBoxFound = cw.subBoxFound;

    if(subBoxFound) {
      SubBox bestSubBox = cw.bestSubBox;
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







