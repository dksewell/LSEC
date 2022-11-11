#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double QUV(const arma::colvec & SS,
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
  int L = IVApprox.nrow();
  
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  double ret = 0;
  
  if(subsamp){
    for(int k=0; k<K; k++){
      for(int i=0;i<n;i++){
        ret += Pmki1(i,k)*(SS(i) + dot(U.row(i),W.row(k))) + 
          Pmki2(i,k)*(RR(i) + dot(V.row(i),W.row(k)));
      }
      ret -= Pk(k)*log(fuk(k));
      for(int ell=0; ell<L; ell++){
        ret -=  Pk(k)*log( fvk(k) - eViWk(EE(IVApprox(ell,k) - 1,0) - 1,k) )/L;
      }
    }
  }else{
    for(int k=0; k<K; k++){
      for(int m = 0;m<M;m++){
        ret += pmk(m,k)*( SS(EE(m,0) - 1) + RR(EE(m,1) - 1) + 
          dot(U.row(EE(m,0) - 1) + V.row(EE(m,1) - 1), W.row(k)) - 
          log(fuk(k)) - log( fvk(k) - eViWk(EE(m,0) - 1,k) ) );
      }
    }
  }
  ret += -0.5*tauU*accu(U%U) - 0.5*tauV*accu(V%V) - 0.5*tauS*sum(SS%SS) - 0.5*tauR*sum(RR%RR);
  
  return ret;
}
