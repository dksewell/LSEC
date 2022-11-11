#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat dQUV(const arma::colvec & SS,
               const arma::colvec & RR,
               const arma::mat & U,
               const arma::mat & V,
               const arma::mat & W,
               const arma::mat & pmk,
               const arma::rowvec & Pk,
               const arma::mat & Pmki1,
               const arma::mat & Pmki2,
               const double & tauS,
               const double & tauR,
               const double & tauU,
               const double & tauV,
               const IntegerMatrix & EE,
               const bool & subsamp,
               const IntegerMatrix & IVApprox){
  int M = EE.nrow();
  int K = W.n_rows;
  int n = U.n_rows;
  int p = U.n_cols;
  int L = IVApprox.nrow();
  
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat gradUV = arma::zeros(n,2*(p + 1));
  arma::mat WTilde = arma::ones(K,p + 1);
  WTilde.cols(1,p) = W;
  
  if(subsamp){
    for(int i=0;i<n;i++){
      for(int k=0; k<K; k++){
        gradUV.row(i).cols(0,p) = gradUV.row(i).cols(0,p) +
          (Pmki1(i,k) - Pk(k)/fuk(k)*eUiWk(i,k))*WTilde.row(k);
        gradUV.row(i).cols(p + 1,2*(p + 1) - 1) = gradUV.row(i).cols(p + 1,2*(p + 1) - 1) +
          Pmki2(i,k)*WTilde.row(k);
        for(int ell=0;ell<L;ell++){
          if(EE(IVApprox(ell,k) - 1,0) - 1 != i){
            gradUV.row(i).cols(p + 1,2*(p + 1) - 1) = gradUV.row(i).cols(p + 1,2*(p + 1) - 1) -
              Pk(k)/L/(fvk(k) - eViWk(EE(IVApprox(ell,k) - 1,0) - 1,k))*eViWk(i,k)*WTilde.row(k);
          }
        }
      }
      
      gradUV(i,0) = gradUV(i,0) - tauS*SS(i);
      gradUV.row(i).cols(1,p) = gradUV.row(i).cols(1,p) - tauU*U.row(i);
      gradUV(i,p + 1) = gradUV(i,p + 1) - tauR*RR(i);
      gradUV.row(i).cols(p + 2,2*(p + 1) - 1) = gradUV.row(i).cols(p + 2,2*(p + 1) - 1) - tauV*V.row(i);
    }
  }else{
    for(int i=0;i<n;i++){
      for(int k=0; k<K; k++){
        gradUV.row(i).cols(0,p) = gradUV.row(i).cols(0,p) +
          (Pmki1(i,k) - Pk(k)/fuk(k)*eUiWk(i,k))*WTilde.row(k);
        gradUV.row(i).cols(p + 1,2*(p + 1) - 1) = gradUV.row(i).cols(p + 1,2*(p + 1) - 1) +
          Pmki2(i,k)*WTilde.row(k);
        for(int m=0;m<M;m++){
          if(EE(m,0) - 1 != i){
            gradUV.row(i).cols(p + 1,2*(p + 1) - 1) = gradUV.row(i).cols(p + 1,2*(p + 1) - 1) -
              Pk(k)/L/(fvk(k) - eViWk(EE(m,0) - 1,k))*eViWk(i,k)*WTilde.row(k);
          }
        }
      }
      
      gradUV(i,0) = gradUV(i,0) - tauS*SS(i);
      gradUV.row(i).cols(1,p) = gradUV.row(i).cols(1,p) - tauU*U.row(i);
      gradUV(i,p + 1) = gradUV(i,p + 1) - tauR*RR(i);
      gradUV.row(i).cols(p + 2,2*(p + 1) - 1) = gradUV.row(i).cols(p + 2,2*(p + 1) - 1) - tauV*V.row(i);
    }
  }
  
  return gradUV;
}
