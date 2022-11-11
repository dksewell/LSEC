#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double QRV(const arma::colvec & RR,
           const arma::mat & V,
           const arma::mat & W,
           const arma::mat & pmk,
           const arma::rowvec & Pk,
           const double & tauR,
           const double & tauV,
           const IntegerMatrix & EE){
  int M = pmk.n_rows;
  int K = W.n_rows;
  
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fvk = sum(eViWk,0);
  double ret = -0.5*tauV*accu(V%V) - 0.5*tauR*sum(RR%RR);
  
  for(int k=0; k<K; k++){
    for(int m=0; m<M; m++){
      ret += pmk(m,k)*( 
        RR(EE(m,1) - 1) + dot(V.row(EE(m,1) - 1), W.row(k)) -
        log(fvk(k) - eViWk(EE(m,0) - 1,k))
      );
    }
  }
  
  return ret;
}
