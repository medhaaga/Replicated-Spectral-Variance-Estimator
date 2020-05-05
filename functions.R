set.seed(1)

####function for creating (x,y,z) coordinates lying on the density plot #######
perspective <- function(A, C, pts)
{
  x = seq(-2, 2*C, length= pts)
  y = x
  f = function(x, y) { 
    r = exp(-.5*(A*(x*y)^2 + x*x + y*y -2*C*x - 2*C*y))
    return(r)
  }
  z = outer(x, y, f)
  z[is.na(z)] = 1
  samples <- list("x"=x, "y"=y, "z"=z)
  return(samples)
}

#######create boomerang MC using Gibbs sampler with params A and C, no. of iterates = nsim, and start######
markov.chain <- function(A,C,nsim,start){
  x.samp <- rep(0,nsim)
  y.samp <- rep(0,nsim)
  
  x.samp[1] <- start[1]
  y.samp[1] <- start[2]
  
  for (i in 2:nsim){
    x.samp[i] <- rnorm(1,mean = C/(1 + A*y.samp[i-1]^2), sd = sqrt(1/(1 + A*y.samp[i-1]^2)))
    y.samp[i] <- rnorm(1,mean = C/(1 + A*x.samp[i]^2), sd = sqrt(1/(1 + A*x.samp[i]^2)))
  }
  chain <- cbind(x.samp, y.samp)
  colnames(chain) <- c("X", "Y")
  return (chain)
}

################find the true mean using adaptIntegrate###############
true.mean <- function(A,C){
  integrand <- function(x) {exp(-.5*(A*(x[1]*x[2])^2 + x[1]*x[1] + x[2]*x[2] -2*C*x[1] - 2*C*x[2]))}
  normalise_const.a <- adaptIntegrate(integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf), tol = 1e-10)
  mean_fn <- function(x) {x[1]*exp(-.5*(A*(x[1]*x[2])^2 + x[1]*x[1] + x[2]*x[2] -2*C*x[1] - 2*C*x[2]))/normalise_const.a$integral}
  expectation <- adaptIntegrate(mean_fn, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf), tol = 1e-10)
  mean.a <- expectation$integral 
  true.mean <- c(mean.a, mean.a)
  return(true.mean)
}

##calculate T^2 statistics for a given sample mean (x), true mean (mean), 
##and estimated covariance matrix - sigma
t2.stat <- function(x,mean,sigma,n){
  return(n*t(mean - x) %*% solve(sigma) %*% (mean - x))
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