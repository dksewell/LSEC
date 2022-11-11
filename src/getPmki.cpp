#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat getPmki(const arma::mat & Pmk,
                  const IntegerVector & Mi1Index,
                  const IntegerMatrix & Mi1){
  int K = Pmk.n_cols;
  int n = Mi1.nrow();
  arma::mat Pmki1 = arma::zeros(n,K);
  
  for(int i=0; i<n; i++){
    if(Mi1Index(i)>0){
      if(Mi1Index(i) == 1){
        Pmki1.row(i) = Pmk.row(Mi1(i,0) - 1);
      }else{
        for(int j=0;j<Mi1Index(i);j++){
          Pmki1.row(i) = Pmki1.row(i) +
            Pmk.row(Mi1(i,j) - 1);
        }
      }
    }
    
  }
  
  return Pmki1;
}
