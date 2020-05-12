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

create.output <- function(A,B,C,m,check.pts,freq,c.prob){
  
  start <- matrix(0, nrow = m, ncol = 2)  #only depends on C
  
  for(i in 1:floor(m/2)){ 
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
  }
  
  critical <- qchisq(c.prob, df=2)
  T.mean <- true.mean(A,B,C)
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    
    asv.samp <- array(0, dim = c(2,2,freq))
    rsv.samp <- array(0, dim = c(2,2,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    for (j in 1:freq){
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      chain <- array(0,dim = c(nsim,2,m))
      sve <- array(0, dim = c(2,2,m))
      rsve <- array(0, dim = c(2,2,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        chain[,,k] <- markov.chain(A,B,C,nsim,start[k,])
        b[k] <- batchSize(chain[,,k], method = "bartlett")
      }
      
      b.avg <- mean(b)
      global.mean <- apply(chain,2,mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
      }
      
      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
      
      chi.sq.asv <- t2.stat(global.mean,T.mean,asv.samp[,,j],nsim*m)
      chi.sq.rsv <- t2.stat(global.mean,T.mean,rsv.samp[,,j],nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp,asv.coverage,rsv.coverage, file = paste(paste("Out/", A, B, C, "/out", m,nsim,A,B,C, sep = "_"),".Rdata", sep = ""))
  }
  
}

##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, A, B, C, m, rep=100){

  start <- matrix(0, nrow = m, ncol = 2)  #only depends on C
  
  for(i in 1:floor(m/2)){ 
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
  }
  
  master.chain.rep <- array(0,dim = c(max, 2, m, rep))
  
  for (r in 1:rep){
    for (k in 1:m){
      master.chain.rep[,,k,r] <- markov.chain(A, B, C, max, start[k,])
    }
  }
  
  conv.pts <- seq(min, max,500)
  l <- length(conv.pts)
  asv.samp <- array(0, dim = c(2,2,l))
  rsv.samp <- array(0, dim = c(2,2,l))
  
  for (j in 1:l){
    nsim = conv.pts[j]
    sve.rep <- array(0, dim = c(2, 2, rep))
    rsve.rep <- array(0, dim = c(2, 2, rep))
    for (r in 1:rep){
      chain <- master.chain.rep[1:nsim,,,r]
      sve <- array(0, dim = c(2,2,m))
      rsve <- array(0, dim = c(2,2,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        b[k] <- batchSize(chain[,,k], method = "bartlett")
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
    }
    
    asv.samp[,,j] <- apply(sve.rep, c(1,2), mean)
    rsv.samp[,,j] <- apply(rsve.rep, c(1,2), mean)
    print(paste("Percentage completion: ", round(100*j/l, 2)))
  }
  save(asv.samp,rsv.samp, file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))
  
}


######################################################################

A <- 2
B <- 9
C <- 7
m = 5

#sims for plotting densities and calculating coverage

check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4)
freq <- 1e3  #100 for now, will change later
rep <- 10
c.prob <- .95
min <- 5e2
max <- 1e5
conv.pts <- seq(min, max, 500)

print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output(A, B, C, m, check.pts, freq, c.prob)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")
convergence(min, max, A, B, C, m, rep)

