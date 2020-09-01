set.seed(1)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
sourceCpp("lag.cpp")
library(mvtnorm)

markov.chain <- function(phi, omega, nsim, start){
  p <- ncol(phi)
  chain <- matrix(0, nrow = nsim, ncol = p)
  chain[1,] <- start

  for (t in 2:nsim){
    chain[t,] <- as.vector(phi%*%chain[t-1,]) + as.vector(rmvnorm(1, mean = rep(0,p), sigma = omega))
  }
  return(chain)
}

target.sigma <- function(phi, omega){
  p <- ncol(phi)
  omega.vec <- omega
  dim(omega.vec) <- c(dim(omega.vec)[1]*dim(omega.vec)[2], 1)
  sigma.vec <- solve(diag(p*p) - kronecker(phi, phi)) %*% omega.vec
  sigma <- matrix(sigma.vec, p, p)
  return(sigma)
}

true.sigma <- function(phi, var){
  p <- ncol(phi)
  sigma <- solve(diag(p) - phi) %*% var + var %*% solve(diag(p) - t(phi)) - var
  return (sigma)
}

##FFT function for calculating RSVe on a centered matrix A
mSVEfft <- function (A, b, method = "bartlett")
{
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag2(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))
  FF <- matrix(0, ncol = p, nrow = 2*n)
  FF[1:n,] <- A
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- mvfftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n, ]) / n )
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- fftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n]) / n )
  }

}


t2.stat <- function(x,mean,sigma,n){
  return(n*t(mean - x) %*% solve(sigma) %*% (mean - x))
}

confidence_interval <- function(vector, interval) {

  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  result <- c("mean" = vec_mean, "lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
