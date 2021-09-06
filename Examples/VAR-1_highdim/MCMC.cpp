#include <RcppArmadillo.h>
#include <math.h>
using namespace arma;
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat markov_chain(arma::mat phi, arma::mat omega, int nsim, arma::vec start)
{
  Function f1("rnorm");
  Function f2("chol");
  Function f3("t");
  int p = phi.n_cols;

  arma::mat a(p,1);
  arma::mat chain(nsim, p);

  // arma::mat omegaSq = f3(f2(Named("x") = omega));
  arma::mat omegaSq = sqrtmat_sympd(omega);
  arma::mat vect(p,1);
  arma::mat b(p,1);
  chain.zeros();
  chain.row(0) = trans(start);

  for (int i=1; i <= (nsim-1); i++)
  {
    a = phi * trans(chain.row(i-1));
    vect = rnorm(p);
    b = omegaSq*vect;
    chain.row(i) =  trans(a+b);
  }
  
return(chain);  
}
