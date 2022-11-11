#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double evalMargPost(const arma::colvec & SS,
                    const arma::colvec & RR,
                    const arma::mat & U,
                    const arma::mat & V,
                    const arma::mat & W,
                    const double & tauS,
                    const double & tauR,
                    const double & tauU,
                    const double & tauV,
                    const arma::rowvec & alph,
                    const IntegerMatrix & EE,
                    const double & a_s,
                    const double & b_s,
                    const double & a_r,
                    const double & b_r,
                    const double & a_u,
                    const double & a_v,
                    const double & b_u,
                    const double & b_v,
                    const double & a_0,
                    const bool & subsamp,
                    const IntegerVector & IVApprox){
  int M = EE.nrow();
  int K = W.n_rows;
  int n = U.n_rows;
  int p = U.n_cols;
  int L = IVApprox.size();
  
  arma::mat eUiWk = exp(SS*arma::ones(1,K) + U*W.t());
  arma::mat eViWk = exp(RR*arma::ones(1,K) + V*W.t());
  arma::rowvec fuk = sum(eUiWk,0);
  arma::rowvec fvk = sum(eViWk,0);
  
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
  
  ret = ret - 0.5*tauS*sum(SS%SS) - 0.5*tauR*sum(RR%RR) - 
    0.5*tauU*accu(U%U) - 0.5*tauV*accu(V%V) - 0.5*accu(W%W) +
    0.5*n*p*(log(tauU) + log(tauV)) +
    (0.5*a_s - 1)*log(tauS) - 0.5*b_s*tauS +
    (0.5*a_r - 1)*log(tauR) - 0.5*b_r*tauR +
    (0.5*a_u - 1)*log(tauU) - 0.5*b_u*tauU +
    (0.5*a_v - 1)*log(tauV) - 0.5*b_v*tauV + (a_0 - 1)*sum(log(alph));
  
  return ret;
}
