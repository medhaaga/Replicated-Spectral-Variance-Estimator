library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
sourceCpp("lag.cpp")


# FFT function for SVE
# Remember that A should be centered
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


n <- 1e2
p <- 2
b <- 10

chain1 <- matrix(rnorm(n*p,mean = 2), ncol = 2, nrow = n)
chain2 <- matrix(rnorm(n*p, mean = 2), ncol = 2, nrow = n)


# ASV estimator
c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
c2.cen.loc <- scale(chain2, scale = FALSE)

sve1 <- mSVEfft(A = c1.cen.loc, b = b)
sve2 <- mSVEfft(A = c2.cen.loc, b = b)

(asv <- (sve1 + sve2)/2)

# RSV estimator
global.mean <- colMeans(rbind(chain1, chain2))

c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
c2.cen <- scale(chain2, center = global.mean, scale =FALSE)

rsve1 <- mSVEfft(A = c1.cen, b = b)
rsve2 <- mSVEfft(A = c2.cen, b = b)

(rsv <- (rsve1 + rsve2)/2)





