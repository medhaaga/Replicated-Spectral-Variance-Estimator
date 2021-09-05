#include <RcppArmadillo.h>
#include <math.h>
using namespace arma;
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix markov_chain(NumericMatrix phi, NumericMatrix omega, int nsim, NumericVector start)
{
  Function f1("rnorm");
  Function f2("chol");
  Function f3("t");
  int p = phi.ncol();
  NumericMatrix chain(nsim, p);
  chain(0,_) = start;
  for (int i=1; i <= (nsim-1); i++){
    NumericVector a = (phi*chain(i-1,_));
    NumericMatrix mat = f3(f2(Named("x") = omega));
    NumericVector vec = f1(p);
    NumericVector b = mat*vec;
    chain(i,_) =  a+b;
  }
  
return(chain);  
}
