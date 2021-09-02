#include <RcppArmadillo.h>
#include <math.h>
using namespace arma;
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix markov_chain(NumericMatrix phi, NumericMatrix omega, int nsim, NumericVector start)
{
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["rmvnorm"];
  int p = phi.ncol();
  NumericMatrix chain(nsim, p);
  chain(0,_) = start;
  for (int i=1; i <= (nsim-1); i++){
    NumericVector a = (phi*chain(i-1,_));
    NumericVector b = f(1, Named("sigma", omega));
    chain(i,_) =  a + b;
  }
  
return(chain);  
}
