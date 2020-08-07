################################################################
##This code is responsible for simulating and storing R objects 
##that are later used for plotting figures and making tables.
################################################################

set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")

############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output <- function(A, B, C, m, check.pts, freq, c.prob){
  
  p <- 2
  start <- matrix(0, nrow = m, ncol = 2)  #only depends on C
  T.mean <- true.mean(A, B, C)
  
  for(i in 1:floor(m/2)){
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
  }
  
  asv.samples <- list()
  gsv.samples <- list()
  asv.coverage <- list()
  gsv.coverage <- list()
  
  for (i in 1:length(check.pts)){
    asv.samples[[i]] <- array(0, dim = c(p,p,freq))
    gsv.samples[[i]] <- array(0, dim = c(p,p,freq))
    asv.coverage[[i]] <- rep(0, freq)
    gsv.coverage[[i]] <- rep(0, freq)
  }
  
  for (j in 1:freq){
    if(j %% (freq/1000) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2)))
    
    master.chain <- array(0, dim = c(max(check.pts), p, m))
    for (k in 1:m)
      master.chain[,,k] <- markov.chain(A, B, C, max(check.pts), start[k,])
    
    for (i in 1:length(check.pts))
    {
      nsim <- check.pts[i]
      critical <- ((nsim*m - 1)*p/(nsim*m  - p))*qf(.95, df1 = p, df2 = (nsim*m-p))
      chains <- master.chain[1:nsim,,]
      
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        b[k] <- batchSize(chains[,,k], method = "bartlett")
      }
      
      b.avg <- ceiling(mean(b))
      global.mean <- apply(chains, 2, mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chains[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chains[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
      }
      asv.samples[[i]][,,j] <- apply(sve, c(1,2), mean)
      gsv.samples[[i]][,,j] <- apply(rsve, c(1,2), mean)
      
      chi.sq.asv <- t2.stat(global.mean, T.mean, asv.samples[[i]][,,j], nsim*m)
      chi.sq.gsv <- t2.stat(global.mean, T.mean, gsv.samples[[i]][,,j], nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[[i]][j]=1}
      if (chi.sq.gsv <= critical) {gsv.coverage[[i]][j]=1}
    }
  }
  save(asv.samples, gsv.samples, asv.coverage, gsv.coverage, 
       file = paste(paste("Out/", A, B, C, "/out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))
  
}


##################################################
#function for storing convergence plots data
##################################################

convergence <- function(min, max, step, A, B, C, m, rep){
  
  p <- 2
  start <- matrix(0, nrow = m, ncol = 2)  #only depends on C
  
  for(i in 1:floor(m/2)){
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
  }
  
  conv.pts <- seq(min, max, step)
  l <- length(conv.pts)
  asv <- list()
  rsv <- list()
  ess.asv <- list()
  ess.rsv <- list()
  
  for (r in 1:rep){
    
    print(paste("Sampling for rep =", r))
    
    asv.samp <- array(0, dim = c(p, p, l))
    rsv.samp <- array(0, dim = c(p, p, l))
    ess.asv.samp <- rep(0, l)
    ess.rsv.samp <- rep(0, l)
    
    master.chain <- array(0, dim = c(max, p, m))
    
    for (k in 1:m){
      print(paste("Sampling chain - ", k))
      master.chain[,,k] <- markov.chain(A, B, C, max, start[k,])
    }
    
    for (j in 1:l){
      
      nsim <- conv.pts[j]
      chain <- master.chain[1:nsim,,]
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      smpl.cov <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        b[k] <- batchSize(chain[,,k], method = "bartlett")
        smpl.cov[,,k] <- cov(chain[,,k])
      }
      
      b.avg <- ceiling(mean(b))
      
      global.mean <- apply(chain,2,mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
        
      }
      
      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
      lambda.rep <- apply(smpl.cov, c(1,2), mean)
      ess.asv.samp[j] <- (det(lambda.rep)/det(asv.samp[,,j]))^(1/p)
      ess.rsv.samp[j] <- (det(lambda.rep)/det(rsv.samp[,,j]))^(1/p)
      
    }
    
    asv[[r]] <- asv.samp
    rsv[[r]] <- rsv.samp
    ess.asv[[r]] <- ess.asv.samp
    ess.rsv[[r]] <- ess.rsv.samp
    
    save(asv, rsv, ess.asv, ess.rsv,
         file = paste(paste("Out/", A, B, C, "/conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))
  }
  
}


#########################################
### Simulation setting-1; A=1, B=3, C=8
#########################################

#### Two chains

A <- 1
B <- 3
C <- 8
m = 2

check.pts <- c(1e3, 5e3, 1e4, 2e4, 5e4, 1e5)
freq <- 1000
rep <- 50
c.prob <- .95
min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)

print("Carrying out 1000 repititions for each value of nsim in check.pts")

t1 <- Sys.time()
create.output(A, B, C, m, check.pts, freq, c.prob)
t2 <- Sys.time()
print(t2-t1)

#### Five chains

m = 5

print("Carrying out 1000 repititions for each value of nsim in check.pts")
t1 <- Sys.time()
create.output(A, B, C, m, check.pts, freq, c.prob)
t2 <- Sys.time()
print(t2-t1)

print("Carrying out simulations for running plots of ASV and RSV in the range(1e3, 1e5")
t1 <- Sys.time()
convergence(min, max, step, A, B, C, m, rep)
t2 <- Sys.time()
print(t2-t1)

#########################################
### Simulation seeting-2; A=1, B=10, C=7
#########################################

#### Two chains

A <- 1
B <- 10
C <- 7
m = 2

check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4, 5e4)
freq <- 1e3  
rep <- 50
c.prob <- .95
min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)

t1 <- Sys.time()
create.output(A, B, C, m, check.pts, freq, c.prob)
t2 <- Sys.time()
print(t2-t1)

#### Five chains


m = 5

print("Carrying out 1000 repititions for each value of nsim in check.pts")
t1 <- Sys.time()
create.output(A, B, C, m, check.pts, freq, c.prob)
t2 <- Sys.time()
print(t2-t1)

print("Carrying out simulations for running plots of ASV and RSV in the range(1e3, 1e5")
t1 <- Sys.time()
convergence(min, max, step, A, B, C, m, rep)
t2 <- Sys.time()
print(t2-t1)
