#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double QW(const arma::colvec & SS,
          const arma::colvec & RR,
          const arma::mat & U,
          const arma::mat & V,
          const arma::mat & W,
          const arma::mat & pmk,
          const arma::rowvec & Pk,
          const IntegerMatrix & EE,
          const bool & subsamp,
          const IntegerMatrix & IVApprox){
  int M = pmk.n_rows;
  int K = W.n_rows;
  int L = IVApprox.nrow();
  
  arma::mat eUiWki = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWki = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWki,0);
  arma::rowvec fvk = sum(eViWki,0);
  double ret = -0.5*accu(W%W);
  
  if(subsamp){
    for(int k=0; k<K; k++){
      for(int ell=0; ell<L; ell++){
        ret += Pk(k)*(
          dot(U.row(EE(IVApprox(ell,k) - 1,0) - 1) + V.row(EE(IVApprox(ell,k) - 1,1) - 1), W.row(k)) -
          log(fuk(k)) - log(fvk(k) - eViWki(EE(IVApprox(ell,k) - 1,0) - 1,k))
        )/L;
      }
    }
  }else{
    for(int k=0; k<K; k++){
      for(int m=0; m<M; m++){
        ret += pmk(m,k)*( 
          dot(U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1), W.row(k)) -
          log(fuk(k)) -
          log(fvk(k) - eViWki(EE(m,0) - 1,k))
        );
      }
    }
  }
  
  return ret;
}

