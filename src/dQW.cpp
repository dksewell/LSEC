#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat dQW(const arma::colvec & SS,
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
  int p = W.n_cols;
  int n = U.n_rows;
  int L = IVApprox.nrow();
  
  arma::mat eUiWki = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWki = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWki,0);
  arma::rowvec fvk = sum(eViWki,0);
  arma::mat gradW = -W;
  
  arma::mat Suk = arma::zeros(K,p);
  arma::mat Svk = arma::zeros(K,p);
  for(int i=0;i<n;i++){
    Suk = Suk + eUiWki.row(i).t()*U.row(i);
    Svk = Svk + eViWki.row(i).t()*V.row(i);
  }
  
  if(subsamp){
    for(int k=0; k<K; k++){
      for(int ell=0; ell<L; ell++){
        gradW.row(k) = gradW.row(k) +
          Pk(k)*(
              U.row(EE(IVApprox(ell,k) - 1,0) - 1) + V.row(EE(IVApprox(ell,k) - 1,1) - 1) -
              Suk.row(k)/fuk(k) - 
              (Svk.row(k) - eViWki(EE(IVApprox(ell,k) - 1,0) - 1,k)*V.row(EE(IVApprox(ell,k) - 1,0) - 1))/
                (fvk(k) - eViWki(EE(IVApprox(ell,k) - 1,0) - 1,k))
          )/L;
      }
    }
  }else{
    for(int k=0; k<K; k++){
      for(int m=0; m<M; m++){
        gradW.row(k) = gradW.row(k) +
          pmk(m,k)*(
              U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1) -
              Suk.row(k)/fuk(k) - 
              (Svk.row(k) - eViWki(EE(m,0) - 1,k)*V.row(EE(m,0) - 1))/
                (fvk(k) - eViWki(EE(m,0) - 1,k))
          );
      }
    }
  }
  
  return gradW;
}

