set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")
library("MCMCpack")
load(file = "mida.rda")


##################################################
#function for storing convergence plots data
##################################################

convergence <- function(min, max, c0, d0, bars, m, rep=100, step=500){
  
  p <- bars+1
  master.chain.rep <- array(0,dim = c(max, p, m, rep))
  
  for (r in 1:rep){
    for (k in 1:m){
      master.chain.rep[,,k,r] <- matrix(MCMCpoissonChange(mida~ 1, m=bars, c0=c0, d0=d0, marginal.likelihood="none", mcmc = max, burnin = 0, thin = 1, seed = round(1e3*runif(1))), nrow = nsim, ncol = p)
    }
    print(paste("No of replications done:", r))
  }
  
  conv.pts <- seq(min, max,step)
  l <- length(conv.pts)
  asv.samp <- array(0, dim = c(p,p,l))
  rsv.samp <- array(0, dim = c(p,p,l))
  ess.asv.samp <- list()
  ess.rsv.samp <- list()
  
  for (j in 1:l){
    nsim = conv.pts[j]
    sve.rep <- array(0, dim = c(p, p, rep))
    rsve.rep <- array(0, dim = c(p, p, rep))
    lambda.rep <- array(0, dim = c(p,p,rep))
    ess.asv.rep <- rep(0,rep)
    ess.rsv.rep <- rep(0,rep)
    
    for (r in 1:rep){
      chain <- master.chain.rep[1:nsim,,,r]
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      smpl.cov <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        b[k] <- batchSize(chain[,,k], method = "bartlett")
        smpl.cov[,,k] <- cov(chain[,,k])
      }
      
      b.avg <- mean(b)
      global.mean <- apply(chain,2,mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
        
      }
      
      sve.rep[,,r] <- apply(sve, c(1,2), mean)
      rsve.rep[,,r] <- apply(rsve, c(1,2), mean)
      lambda.rep[,,r] <- apply(smpl.cov, c(1,2), mean)
      ess.asv.rep[r] <- (det(lambda.rep[,,r])/det(sve.rep[,,r]))^(1/p)
      ess.rsv.rep[r] <- (det(lambda.rep[,,r])/det(rsve.rep[,,r]))^(1/p)
    }
    
    asv.samp[,,j] <- apply(sve.rep, c(1,2), mean)
    rsv.samp[,,j] <- apply(rsve.rep, c(1,2), mean)
    ess.asv.samp[[j]] <- ess.asv.rep
    ess.rsv.samp[[j]] <- ess.rsv.rep
    print(paste("Percentage completion: ", round(100*j/l, 2)))
  }
  save(asv.samp,rsv.samp, ess.asv.samp, ess.rsv.samp, file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))
  
}


######################################
###########Model Parameters############
######################################

c0 <- 13
d0 <- 1
bars <- 6
m <- 2       #### Number of chains

#sims for plotting run plots

p <- bars+1  #### dimension
rep <- 10
min <- 5e2
max <- 1e5
step <- 500
conv.pts <- seq(min, max, step)
load(file = "true_mida_1e6_1e3")
truth <- colMeans(mida_means)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")
convergence(min, max, c0, d0, bars, m, rep, step=500)

