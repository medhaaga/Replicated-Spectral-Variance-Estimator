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
  
  start <- matrix(c(2*C,1,1,2*C), 2, 2)  #only depends on C
  critical <- qchisq(c.prob, df=2)
  T.mean <- true.mean(A,B,C)
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    print(paste("Simulationg for nsim = ", nsim))
    
    asv.samp <- array(0, dim = c(2,2,freq))
    rsv.samp <- array(0, dim = c(2,2,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    for (j in 1:freq){
      print(paste("Percentage completion: ", round(j/10, 2)))
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
    save(asv.samp,rsv.samp,asv.coverage,rsv.coverage, file = paste(paste("out",nsim,A,B,C, sep = "_"),".Rdata", sep = ""))
  }
  
}

##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, A, B, C){

  start <- matrix(c(2*C,1,1,2*C), 2, 2)
  master.chain <- array(0,dim = c(max,2,m))
  
  for (k in 1:m){
    master.chain[,,k] <- markov.chain(A,B,C,max,start[k,])
  }
  conv.pts <- seq(min, max,500)
  
  l <- length(conv.pts)
  asv.samp <- array(0, dim = c(2,2,l))
  rsv.samp <- array(0, dim = c(2,2,l))
  
  for (j in 1:l){
    nsim = conv.pts[j]
    chain <- master.chain[1:nsim,,]
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
    
    asv.samp[,,j] <- apply(sve, c(1,2), mean)
    rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
    print(paste("Percentage completion: ", round(100*j/l, 2)))
  }
  save(asv.samp,rsv.samp, file = paste(paste("conv_data", min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))
  
}


######################################################################
#10*3 matrix for A,B,C parameters. Each row corresponds to a different suitable parameterization
params <- matrix(c(1,2,7,1,8,9,1,9,9,2,6,7,2,8,7,2,9,7,2,10,8), nrow = 7, ncol=3, byrow = TRUE)
m = 2
#sims for plotting densities and calculating coverage

check.pts <- c(1e3, 5e3, 1e4, 5e4, 7e4, 1e5)
freq <- 1e3
c.prob <- .95
min <- 1e3
max <- 1e5

# for creating .Rdata files for each set of parameter values
for (t in 1:7)
{
  print(paste("Sampling for A, B, C, = ", params[t,1], params[t,2], params[t,3], "respectively", sep = " "))
  print("Carrying out 1000 repititions for each value of nsim in check.pts")
  create.output(params[t,1], params[t,2], params[t,3], m, check.pts, freq, c.prob)
  print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")
  convergence(min, max, params[t,1], params[t,2], params[t,3])
}
