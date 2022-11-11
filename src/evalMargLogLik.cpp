#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double evalMargLogLik(const arma::colvec & SS,
                      const arma::colvec & RR,
                      const arma::mat & U,
                      const arma::mat & V,
                      const arma::mat & W,
                      const arma::rowvec & alph,
                      const IntegerMatrix & EE,
                      const bool & subsamp,
                      const IntegerVector & IVApprox){
  int M = EE.nrow();
  int K = W.n_rows;
  int L = IVApprox.size();
  
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
  arma::mat Pmk = arma::zeros(M,K);
  arma::colvec PmkRowSums = arma::zeros(M);
  
  double ret = 0;
  
  if(subsamp){
    arma::mat Pmk = arma::zeros(L,K);
    arma::colvec PmkRowSums = arma::zeros(L);
    
    for(int ell = 0;ell<L;ell++){
      Pmk.row(ell) = 
        eUiWk.row(EE(IVApprox(ell) - 1,0) - 1)%eViWk.row(EE(IVApprox(ell) - 1,1) - 1)%alph/
          fuk/(fvk - eViWk.row(EE(IVApprox(ell) - 1,0) - 1));
    }
    PmkRowSums = log(sum(Pmk,1));
    ret = sum(PmkRowSums)/L*M;
  }else{
    arma::mat Pmk = arma::zeros(M,K);
    arma::colvec PmkRowSums = arma::zeros(M);
    
    for(int m = 0;m<M;m++){
      Pmk.row(m) = 
        eUiWk.row(EE(m,0) - 1)%eViWk.row(EE(m,1) - 1)%alph/fuk/(fvk - eViWk.row(EE(m,0) - 1));
    }
    PmkRowSums = log(sum(Pmk,1));
    ret = sum(PmkRowSums);
  }
  
  return ret;
}

