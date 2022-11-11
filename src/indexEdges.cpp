#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List indexEdges(const IntegerMatrix & EE,
                      const int & n){
  IntegerMatrix Mi1(n,n);
  IntegerMatrix Mi2(n,n);
  IntegerVector Mi1Index(n);
  IntegerVector Mi2Index(n);
  int M = EE.nrow();
  
  for(int m=0;m<M;m++){
    Mi1(EE(m,0) - 1, Mi1Index(EE(m,0) - 1)) = m + 1;
    Mi2(EE(m,1) - 1, Mi2Index(EE(m,1) - 1)) = m + 1;
    Mi1Index(EE(m,0) - 1) = Mi1Index(EE(m,0) - 1) + 1;
    Mi2Index(EE(m,1) - 1) = Mi2Index(EE(m,1) - 1) + 1;
  }
  
  return List::create(Named("Mi1") = Mi1,
                      Named("Mi2") = Mi2,
                      Named("Mi1Index") = Mi1Index,
                      Named("Mi2Index") = Mi2Index);
}

