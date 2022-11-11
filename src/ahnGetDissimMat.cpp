#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat ahnGetDissimMat(const arma::mat & simMatij,
                          const IntegerMatrix & EL){
  int M = EL.nrow();
  arma::mat dMat = arma::zeros(M,M);
  
  /*
  for(int m1=0;m1<M - 1;m1++){
  for(int m2=m1 + 1;m2<M;m2++){
  if(EL(m1,1) == EL(m2,1)){
  dMat(m1,m2) = dMat(m2,m1) = 
  1.0 - simMatij(EL(m1,0) - 1, EL(m2,0) - 1);
  }
  }
  }
  */
  
  for(int m1=0;m1<M - 1;m1++){
    for(int m2=m1 + 1;m2<M;m2++){
      if( (EL(m1,0) == EL(m2,0)) ){
        dMat(m1,m2) = dMat(m2,m1) = 
          1.0 - simMatij(EL(m1,1) - 1, EL(m2,1) - 1);
      }
      if( (EL(m1,0) == EL(m2,1)) ){
        dMat(m1,m2) = dMat(m2,m1) = 
          1.0 - simMatij(EL(m1,1) - 1, EL(m2,0) - 1);
      }
      if( (EL(m1,1) == EL(m2,0)) ){
        dMat(m1,m2) = dMat(m2,m1) = 
          1.0 - simMatij(EL(m1,0) - 1, EL(m2,1) - 1);
      }
      if( (EL(m1,1) == EL(m2,1)) ){
        dMat(m1,m2) = dMat(m2,m1) = 
          1.0 - simMatij(EL(m1,0) - 1, EL(m2,0) - 1);
      }
    }
  }
  
  return dMat;
}

