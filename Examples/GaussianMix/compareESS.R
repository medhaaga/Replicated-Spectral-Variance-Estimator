set.seed(1)
library(multichainACF)


log.density <- function(x, p, mu1, mu2, sd1, sd2){
  return(log(p*dnorm(x, mean = mu1, sd = sd1) + (1-p)*dnorm(x, mean = mu2, sd= sd2)))
}


mh.mcmc <- function(start, p, mu1, mu2, sd1, sd2, N, h){
  X <- rep(0, N)
  X[1] <- start
  for (i in 2:N){
    prop <- rnorm(1, mean = X[i-1], sd = 2)
    ratio <- log.density(prop, p, mu1, mu2, sd1, sd2) - log.density(X[i-1], p, mu1, mu2, sd1, sd2)
    if(runif(1) < exp(ratio))
      X[i] <- prop else
        X[i] <- X[i-1]
  }
  return (X)
}


N <- 1e5
p <- 0.7
mu1 <- -5
mu2 <- 5
sd1 <- 1
sd2 <- 0.5


chain1 <- mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N, 1)
chain2 <- mh.mcmc(start = 3, p, mu1, mu2, sd1, sd2, N, 1)
mean(chain1)
mean(chain2)


library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")



ind <- seq(500, N, by = 500)

ess1 <- numeric(length = length(ind))
ess2 <- numeric(length = length(ind))
for(i in 1:length(ind))
{
  if(i%%10 == 0) print(i/length(ind)*100)
  chain1.sub <- chain1[1:ind[i]]
  chain2.sub <- chain2[1:ind[i]]

  chains <- cbind(chain1.sub, chain2.sub)
  b.avg <- mean(batchSize(chain1.sub, method = "bartlett"), batchSize(chain2.sub, method = "bartlett") )

  global.mean <- mean(chains)
  chain.cen.loc <- scale(chains, scale = FALSE)  ## X_st - bar(X)_s
  sve <- mean(diag(mSVEfft(A = chain.cen.loc, b = b.avg)))
  chain.cen <- chains - global.mean
  rsve <- mean(diag(mSVEfft(A = chain.cen, b = b.avg)))



  ess1[i] <- var(c(chain1.sub,chain2.sub))/rsve
  ess2[i] <- mean(apply(chains, 2, var))/rsve  
}

pdf("estimated_ESS.pdf", height = 4, width = 5)
plot(ind, ess1, ylim = range(c(ess1,ess2)),type = 'l', col = "red", ylab = "Estimated ESS/n", xlab = "Sample Size", cex.lab = 1.2, cex.axis = 1.2)
lines(ind, ess2, type = 'l', col = "blue")
dev.off()

