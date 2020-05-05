set.seed(1)

library(fields)
library(scatterplot3d)
library(cubature)
library(graphics)
library(pracma)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
sourceCpp("lag.cpp")
source_python("num_integration.py")  #comment it not using a python environment
scipy <- import("scipy")

###############################################
#####Visualization of distribution############
###############################################

load(file = "perspective.Rdata")
#Perspective plot
persp(x,y,z,col="lightblue",main="Perspective Plot")

#Contour plot
pdf(file = "contour_plot.pdf")
contour(x,y,z,main="Contour Plot")
filled.contour(x,y,z,color=terrain.colors,main="Contour Plot",)
dev.off()

######################################
###Finding the normalising constant and expectation values###
######################################

# a.) Using num integration from pracma:integral
integrand <- function(x) {exp(-.5*(A*(x[1]*x[2])^2 + x[1]*x[1] + x[2]*x[2] -2*C*x[1] - 2*C*x[2]))}
normalise_const.a <- adaptIntegrate(integrand, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))
mean_fn <- function(x) {x[1]*exp(-.5*(A*(x[1]*x[2])^2 + x[1]*x[1] + x[2]*x[2] -2*C*x[1] - 2*C*x[2]))/normalise_const.a$integral}
expectation <- adaptIntegrate(mean_fn, lowerLimit = c(-Inf, -Inf), upperLimit = c(Inf, Inf))
mean.a <- expectation$integral 
true.mean.r <- c(mean.a, mean.a)
 
# b.) Using numerical integration from scipy. Comment if not using a python environment
normalise_const.b <- num_integration(A,C)
mean.b <- expectation_x(A,C,normalise_const.b[1])[1]
true.mean.py <- c(as.numeric(mean.b), as.numeric(mean.b))

###############################################
###FFT SVE calculation function
###############################################

#A is centered
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

################################################
#########Comaprison of RSV and ASV#############
###############################################

# For n = 1000

n1 <- 1000

x1.samp <- rep(0,n1)
y1.samp <- rep(0,n1)
x2.samp <- rep(0,n1)
y2.samp <- rep(0,n1)

x1.samp[1] <- 15
y1.samp[1] <- 1
x2.samp[1] <- 1
y2.samp[1] <- 15

for (i in 2:n1){
  x1.samp[i] <- rnorm(1,mean = C/(1 + y1.samp[i-1]^2), sd = sqrt(1/(1 + y1.samp[i-1]^2)))
  y1.samp[i] <- rnorm(1,mean = C/(1 + x1.samp[i]^2), sd = sqrt(1/(1 + x1.samp[i]^2)))
  x2.samp[i] <- rnorm(1,mean = C/(1 + y2.samp[i-1]^2), sd = sqrt(1/(1 + y2.samp[i-1]^2)))
  y2.samp[i] <- rnorm(1,mean = C/(1 + x2.samp[i]^2), sd = sqrt(1/(1 + x2.samp[i]^2)))
}

pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Markov-chains-scatter-plot_n1.pdf")
plot(x1.samp, y1.samp, col = rep("red",n1), xlim = c(-2,20), ylim = c(-2,20), main = paste("Scatter plot of two Markov chains, n = ", n1), xlab = "X", ylab = "Y")
points(x2.samp, y2.samp, col = rep("green",n1))
dev.off()

chain1 <- cbind(x1.samp, y1.samp)
chain2 <- cbind(x2.samp, y2.samp)
colnames(chain1) <- c("X", "Y")
colnames(chain2) <- c("X", "Y")

#truncation point 
b1 <- batchSize(chain1, method = "bartlett")
b2 <- batchSize(chain2, method = "bartlett")
b <- (b1+b2)/2

#chain centered by their individual mean
c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
c2.cen.loc <- scale(chain2, scale = FALSE)

#chains centered by global mean
global.mean.n1 <- colMeans(rbind(chain1, chain2))
c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
c2.cen <- scale(chain2, center = global.mean, scale =FALSE)


#ASV estimator

sve1 <- mSVEfft(A = c1.cen.loc, b = b)
sve2 <- mSVEfft(A = c2.cen.loc, b = b)

(asv_n1 <- (sve1 + sve2)/2)

# RSV estimator

rsve1 <- mSVEfft(A = c1.cen, b = b)
rsve2 <- mSVEfft(A = c2.cen, b = b)

(rsv_n1 <- (rsve1 + rsve2)/2)

N <- 100
rsv.samp <- array(0, dim = c(2,2,N))
asv.samp <- array(0, dim = c(2,2,N))

for (j in 1:N){
  x1.samp[1] <- 15
  y1.samp[1] <- 1
  x2.samp[1] <- 1
  y2.samp[1] <- 15
  
  for (i in 2:n1){
    x1.samp[i] <- rnorm(1,mean = C/(1 + y1.samp[i-1]^2), sd = sqrt(1/(1 + y1.samp[i-1]^2)))
    y1.samp[i] <- rnorm(1,mean = C/(1 + x1.samp[i]^2), sd = sqrt(1/(1 + x1.samp[i]^2)))
    x2.samp[i] <- rnorm(1,mean = C/(1 + y2.samp[i-1]^2), sd = sqrt(1/(1 + y2.samp[i-1]^2)))
    y2.samp[i] <- rnorm(1,mean = C/(1 + x2.samp[i]^2), sd = sqrt(1/(1 + x2.samp[i]^2)))
  }
  chain1 <- cbind(x1.samp, y1.samp)
  chain2 <- cbind(x2.samp, y2.samp)
  b1 <- batchSize(chain1, method = "bartlett")
  b2 <- batchSize(chain2, method = "bartlett")
  b <- (b1+b2)/2

  c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
  c2.cen.loc <- scale(chain2, scale = FALSE)

  global.mean <- colMeans(rbind(chain1, chain2))
  c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
  c2.cen <- scale(chain2, center = global.mean, scale =FALSE)
  sve1 <- mSVEfft(A = c1.cen.loc, b = b)
  sve2 <- mSVEfft(A = c2.cen.loc, b = b)
  asv.samp[,,j] <- (sve1 + sve2)/2
  
  rsve1 <- mSVEfft(A = c1.cen, b = b)
  rsve2 <- mSVEfft(A = c2.cen, b = b)
  rsv.samp[,,j] <- (rsve1+rsve2)/2
}

rsv.sigma11 <- rsv.samp[1,1,]
asv.sigma11 <- asv.samp[1,1,]
density.rsv.sigma11 <- stats::density(rsv.sigma11)
density.asv.sigma11 <- stats::density(asv.sigma11)
rsv.sigma22 <- rsv.samp[2,2,]
asv.sigma22 <- asv.samp[2,2,]
density.rsv.sigma22 <- stats::density(rsv.sigma22)
density.asv.sigma22 <- stats::density(asv.sigma22)

####################Variance of Y-component for RSV and ASV#############################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Sigma_11_n1.pdf")
plot(density.rsv.sigma11, col = "red", xlab = "Variance of X-component", main = paste("Density plot for variation of X-component, N=",N))
lines(density.asv.sigma11, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
########################################################################################
####################Variance of Y-component for RSV and ASV#############################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Sigma_22_n1.pdf")
plot(density.rsv.sigma22, col = "red", xlab = "Variance of Y-component", main = paste("Density plot for variation of Y-component, N=",N))
lines(density.asv.sigma22, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
#########################################################################################

det.rsv <- rep(0,N)
det.asv <- rep(0,N)
for (i in 1:N){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}
density.rsv.det <- stats::density(det.rsv)
density.asv.det <- stats::density(det.asv)

####################determinant comparison for RSV and ASV################################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/determinant_n1.pdf")
plot(density.rsv.det, col = "red", xlab = "Determinant of covariance matrix", main = paste("Density plot for determinant of estimated variance, N=",N))
lines(density.asv.sigma22, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
##########################################################################################


#n=50000. We hope to that the Markov chiains explore both modes here

n2 <- 50000

x1.samp <- rep(0,n2)
y1.samp <- rep(0,n2)
x2.samp <- rep(0,n2)
y2.samp <- rep(0,n2)

x1.samp[1] <- 15
y1.samp[1] <- 1
x2.samp[1] <- 1
y2.samp[1] <- 15

for (i in 2:n2){
  x1.samp[i] <- rnorm(1,mean = C/(1 + y1.samp[i-1]^2), sd = sqrt(1/(1 + y1.samp[i-1]^2)))
  y1.samp[i] <- rnorm(1,mean = C/(1 + x1.samp[i]^2), sd = sqrt(1/(1 + x1.samp[i]^2)))
  x2.samp[i] <- rnorm(1,mean = C/(1 + y2.samp[i-1]^2), sd = sqrt(1/(1 + y2.samp[i-1]^2)))
  y2.samp[i] <- rnorm(1,mean = C/(1 + x2.samp[i]^2), sd = sqrt(1/(1 + x2.samp[i]^2)))
}

pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Markov-chains-scatter-plot_n2.pdf")
plot(x1.samp, y1.samp, col = rep("red",n2), xlim = c(-2,20), ylim = c(-2,20), main = paste("Scatter plot of two Markov chains, n = ", n2), xlab = "X", ylab = "Y")
points(x2.samp, y2.samp, col = rep("green",n2))
dev.off()

chain1 <- cbind(x1.samp, y1.samp)
chain2 <- cbind(x2.samp, y2.samp)
colnames(chain1) <- c("X", "Y")
colnames(chain2) <- c("X", "Y")

#truncation point 
b1 <- batchSize(chain1, method = "bartlett")
b2 <- batchSize(chain2, method = "bartlett")
b <- (b1+b2)/2

#chain centered by their individual mean
c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
c2.cen.loc <- scale(chain2, scale = FALSE)

#chains centered by global mean
global.mean.n2 <- colMeans(rbind(chain1, chain2))
c1.cen <- scale(chain1, center = global.mean.n2, scale =FALSE)
c2.cen <- scale(chain2, center = global.mean.n2, scale =FALSE)


#ASV estimator

sve1 <- mSVEfft(A = c1.cen.loc, b = b)
sve2 <- mSVEfft(A = c2.cen.loc, b = b)

(asv_n2 <- (sve1 + sve2)/2)

# RSV estimator

rsve1 <- mSVEfft(A = c1.cen, b = b)
rsve2 <- mSVEfft(A = c2.cen, b = b)

(rsv_n2 <- (rsve1 + rsve2)/2)

N <- 100
rsv.samp <- array(0, dim = c(2,2,N))
asv.samp <- array(0, dim = c(2,2,N))

for (j in 1:N){
  x1.samp[1] <- 15
  y1.samp[1] <- 1
  x2.samp[1] <- 1
  y2.samp[1] <- 15
  
  for (i in 2:n2){
    x1.samp[i] <- rnorm(1,mean = C/(1 + y1.samp[i-1]^2), sd = sqrt(1/(1 + y1.samp[i-1]^2)))
    y1.samp[i] <- rnorm(1,mean = C/(1 + x1.samp[i]^2), sd = sqrt(1/(1 + x1.samp[i]^2)))
    x2.samp[i] <- rnorm(1,mean = C/(1 + y2.samp[i-1]^2), sd = sqrt(1/(1 + y2.samp[i-1]^2)))
    y2.samp[i] <- rnorm(1,mean = C/(1 + x2.samp[i]^2), sd = sqrt(1/(1 + x2.samp[i]^2)))
  }
  chain1 <- cbind(x1.samp, y1.samp)
  chain2 <- cbind(x2.samp, y2.samp)
  b1 <- batchSize(chain1, method = "bartlett")
  b2 <- batchSize(chain2, method = "bartlett")
  b <- (b1+b2)/2
  
  c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
  c2.cen.loc <- scale(chain2, scale = FALSE)
  
  global.mean <- colMeans(rbind(chain1, chain2))
  c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
  c2.cen <- scale(chain2, center = global.mean, scale =FALSE)
  sve1 <- mSVEfft(A = c1.cen.loc, b = b)
  sve2 <- mSVEfft(A = c2.cen.loc, b = b)
  asv.samp[,,j] <- (sve1 + sve2)/2
  
  rsve1 <- mSVEfft(A = c1.cen, b = b)
  rsve2 <- mSVEfft(A = c2.cen, b = b)
  rsv.samp[,,j] <- (rsve1+rsve2)/2
}

rsv.sigma11 <- rsv.samp[1,1,]
asv.sigma11 <- asv.samp[1,1,]
density.rsv.sigma11 <- stats::density(rsv.sigma11)
density.asv.sigma11 <- stats::density(asv.sigma11)
rsv.sigma22 <- rsv.samp[2,2,]
asv.sigma22 <- asv.samp[2,2,]
density.rsv.sigma22 <- stats::density(rsv.sigma22)
density.asv.sigma22 <- stats::density(asv.sigma22)

####################Variance of Y-component for RSV and ASV#############################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Sigma_11_n2.pdf")
plot(density.rsv.sigma11, col = "red", xlab = "Variance of X-component", main = paste("Density plot for variation of X-component, N=",N))
lines(density.asv.sigma11, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
########################################################################################
####################Variance of Y-component for RSV and ASV#############################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/Sigma_22_n2.pdf")
plot(density.rsv.sigma22, col = "red", xlab = "Variance of Y-component", main = paste("Density plot for variation of Y-component, N=",N))
lines(density.asv.sigma22, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
#########################################################################################

det.rsv <- rep(0,N)
det.asv <- rep(0,N)
for (i in 1:N){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}
density.rsv.det <- stats::density(det.rsv)
density.asv.det <- stats::density(det.asv)

####################determinant comparison for RSV and ASV################################
pdf(file = "C:/Users/DELL/Documents/GitHub/Replicated-Spectral-Variance-Estimator/determinant_n2.pdf")
plot(density.rsv.det, col = "red", xlab = "Determinant of covariance matrix", main = paste("Density plot for determinant of estimated variance, N=",N))
lines(density.asv.sigma22, col = "blue")
legend("topright", legend=c("RSV", "ASV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()
##########################################################################################

###########################################
##Coverage Proababilities for sample size n*m
############################################

hotelling_t2_stat <- function(x,mean,sigma,n){
  return(n*t(mean - x) %*% solve(sigma) %*% (mean - x))
}
critical.95.n1 <- qf(.95, df1 = 2, df2 = n1-2)
coverage.rsv <- 0
coverage.asv <- 0
n <- 1000

x1.samp <- rep(0,n1)
y1.samp <- rep(0,n1)
x2.samp <- rep(0,n1)
y2.samp <- rep(0,n1)

for (j in 1:n){
  x1.samp[1] <- 15
  y1.samp[1] <- 1
  x2.samp[1] <- 1
  y2.samp[1] <- 15
  
  for (i in 2:n1){
    x1.samp[i] <- rnorm(1,mean = C/(1 + y1.samp[i-1]^2), sd = sqrt(1/(1 + y1.samp[i-1]^2)))
    y1.samp[i] <- rnorm(1,mean = C/(1 + x1.samp[i]^2), sd = sqrt(1/(1 + x1.samp[i]^2)))
    x2.samp[i] <- rnorm(1,mean = C/(1 + y2.samp[i-1]^2), sd = sqrt(1/(1 + y2.samp[i-1]^2)))
    y2.samp[i] <- rnorm(1,mean = C/(1 + x2.samp[i]^2), sd = sqrt(1/(1 + x2.samp[i]^2)))
  }
  chain1 <- cbind(x1.samp, y1.samp)
  chain2 <- cbind(x2.samp, y2.samp)
  mean.statistic <- colMeans(rbind(chain1, chain2))
  t2_stat.rsv <- hotelling_t2_stat(mean.statistic,true.mean.py,rsv_n1,2*n1)   ##alternatively use true mean 
  t2_stat.asv <- hotelling_t2_stat(mean.statistic,true.mean.py,asv_n1,2*n1)   ##calculated using pracma
  
  if (t2_stat.rsv <= critical.95.n1) {
    coverage.rsv = coverage.rsv+1
    }
  if (t2_stat.asv <= critical.95.n1){ 
    coverage.asv = coverage.asv+1
    }
}
(coverage.rsv/n)
(coverage.asv/n)

