library(timedelay)

markov.chain <- function(nsim, delta.start, lcA, lcB, micro){

  chain <- bayesian(data.lcA = lcA, data.lcB = lcB,
                           data.flux = FALSE, theta.ini = c(0, 0.01, 200),
                           delta.ini = delta.start, delta.uniform.range = c(0, 100),
                           delta.proposal.scale = 1,
                           tau.proposal.scale = 3,
                           tau.prior.shape = 1, tau.prior.scale = 1,
                           sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7,
                           asis = TRUE, micro = micro,
                           sample.size = nsim, warmingup.size = 0, adaptive.delta = FALSE,
                           adaptive.tau = FALSE)
  samples <- cbind(chain$delta, chain$beta, chain$X, chain$mu,
                        chain$sigma, chain$tau)
  chain <- list("samples" = samples, "delta.accept.rate" = chain$delta.accept.rate,
                       "tau.accept.rate" = chain$tau.accept.rate)
  return(chain)
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

confidence_interval <- function(vector, interval) {
  
  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  result <- c("mean" = vec_mean, "lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}
