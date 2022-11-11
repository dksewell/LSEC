#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat computePmk(const arma::colvec & SS,
                     const arma::colvec & RR,
                     const arma::mat & U,
                     const arma::mat & V,
                     const arma::mat & W,
                     const arma::rowvec & alph,
                     const IntegerMatrix & EE){
  int M = EE.nrow();
  int K = W.n_rows;
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat Pmk = arma::zeros(M,K);
  
  for(int m = 0;m<M;m++){
    Pmk.row(m) = 
      eUiWk.row(EE(m,0) - 1)%eViWk.row(EE(m,1) - 1)%alph/fuk/(fvk - eViWk.row(EE(m,0) - 1));
    Pmk.row(m) = Pmk.row(m)/sum(Pmk.row(m));
  }
  
  return Pmk;
}
