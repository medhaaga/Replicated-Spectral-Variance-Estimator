set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
sourceCpp("lag.cpp")
library(spBayes)
data("NETemp.dat")
ne.temp <- NETemp.dat


spatio_temp <- function(N = 5e4, verbose = FALSE, starting)
{
  ##take a chunk of New England
  ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]##subset first 2 years (Jan 2000 - Dec. 2002)
  y.t <- ne.temp[,4:15]
  N.t <- ncol(y.t) ##number of months
  n <- nrow(y.t) ##number of observation per months
  
  coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
  max.d <- max(iDist(coords))
  
  ##set starting and priors
  p <- 2 #number of regression parameters in each month
  
  tuning <- list("phi"=rep(3, N.t))
  
  
  priors <- list("beta.0.Norm"=list(rep(0,p), diag(400,p)),
                 "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
                 "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
                 "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
                 "sigma.eta.IW"=list(2, diag(0.001,p)))
  
  
  ##make symbolic model formula statement for each month
  mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)
  n.samples <- N
  
  m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
                 starting=starting, tuning=tuning, priors=priors, get.fitted =FALSE,
                 cov.model="exponential", n.samples=n.samples, n.report=1e4, verbose = verbose)
  
  
  chain <- cbind(m.1$p.beta.0.samples,
                 m.1$p.beta.samples,
                 m.1$p.theta.samples,
                 m.1$p.sigma.eta.samples[ ,-3],
                 t(m.1$p.u.samples))
  
  
  foo <- rep(1:N.t, each = dim(y.t)[1])
  foo <- paste("t",foo, "s", rep(1:dim(y.t)[1], N.t), sep = "")
  foo <- paste("u.", foo, sep = "")
  colnames(chain)[(p + N.t*p + 3*N.t + p*(p+1)/2 + 1):dim(chain)[2]] <- foo
  
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

