#include <Rcpp.h>
using namespace Rcpp;

//
// Smoothing of distance matrix
//
// param dmat distance matrix (objects in rows and columns)
// param umat membership matrix (objects in rows, clusters in columns)
//
// [[Rcpp::export(".distanceMatrixSmoothing")]]
NumericMatrix distanceMatrixSmoothing(NumericMatrix dmat, NumericMatrix umat, bool add = TRUE) {
  int N = dmat.nrow();
  int K = umat.ncol();
  NumericVector card(K, 0.0);
  NumericVector var(K, 0.0);
  for(int k=0;k<K;k++){
    card[k] = sum(umat(_,k));
    double vg = 0.0;
    for(int i1=0;i1<N;i1++) {
      for(int i2=0;i2<N;i2++) {
        vg += umat(i1,k)*umat(i2,k)*pow(dmat(i1,i2),2.0);
      }
    }
    var[k] = vg/(2.0*pow(card[k], 2.0));
  }
  NumericMatrix dbc(K,K);
  for(int k1=0;k1<K;k1++){
    dbc(k1, k1) = 0.0;
    for(int k2=0;k2<k1;k2++){
      double cd = 0.0;
      for(int i1=0;i1<N;i1++) {
        for(int i2=0;i2<N;i2++) {
          cd += umat(i1,k1)*pow(dmat(i1,i2),2.0)*umat(i2,k2);
        }
      }
      double dsquared = (cd/(card[k1]*card[k2])) - var[k1] - var[k2];
      if(dsquared >= 0.0) {
        dbc(k1,k2) = sqrt(dsquared);
      } else {
        if(add) dbc(k1,k2) = 0.0;
        else dbc(k1,k2) = NA_REAL;
      }
      dbc(k2,k1) = dbc(k1,k2);
    }
  }
  return(dbc);
}

// Smoothing of a rectangular matrix (i.e. multivariate averages)
//
// param smat rectangular data matrix (objects in rows and variables in columns)
// param umat membership matrix (objects in rows, clusters in columns)
// [[Rcpp::export(".rectangularMatrixSmoothing")]]
NumericMatrix rectangularMatrixSmoothing(NumericMatrix smat, NumericMatrix umat) {
  int N = smat.nrow();
  int P = smat.ncol();
  if(umat.nrow()!=umat.ncol()) stop("umat should be a squared matrix");
  if(umat.ncol()!=N) stop("umat should have the same number of rows as smat");
  NumericVector card(N, 0.0);
  for(int i=0;i<N;i++){
    card[i] = sum(umat(_,i));
  }
  NumericMatrix smatout(N,P);
  smatout.fill(0.0);
  for(int i=0;i<N;i++) {
    for(int k=0;k<N;k++) {
      smatout(i,_) = smatout(i,_) + smat(k,_)*umat(k,i);
    }
    smatout(i,_) = smatout(i,_)/card[i];
  }
  return(smatout);
}
