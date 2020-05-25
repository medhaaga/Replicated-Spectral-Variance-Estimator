set.seed(10)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
library(sandwich)
sourceCpp("lag.cpp")

############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output1 <- function(phi, omega, m, check.pts, freq, c.prob){

  p <- ncol(phi)
  truth <- target.sigma(phi, omega)
  start <- matrix(0, nrow = m, ncol = p)  #only depends on C

  for(i in floor(m/2):1){
    start[i,] <- 2*i*sqrt(diag(truth))
    start[m-i+1,] <- -2*i*sqrt(diag(truth))
  }

  critical <- qchisq(c.prob, df=2)

  for (i in 1:length(check.pts)){

    nsim <- check.pts[i]

    asv.samp <- array(0, dim = c(p,p,freq))
    rsv.samp <- array(0, dim = c(p,p,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)

    for (j in 1:freq){
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      chain <- array(0,dim = c(nsim,p,m))
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)

      for (k in 1:m){
        chain[,,k] <- markov.chain(phi, omega, nsim, start[k,])
        b[k] <- batchSize(chain[,,k], method = "bartlett")
      }

      b.avg <- mean(b)
      global.mean <- apply(chain, 2, mean)

      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
      }

      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)

      chi.sq.asv <- t2.stat(global.mean, rep(0,p), asv.samp[,,j], nsim*m)
      chi.sq.rsv <- t2.stat(global.mean, rep(0,p), rsv.samp[,,j], nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp,asv.coverage,rsv.coverage, file = paste(paste("Out/bartlett/out", m, nsim, sep = "_"),".Rdata", sep = ""))
  }

}

##################################################
#function for storing convergence plots data
##################################################
convergence1 <- function(min, max, phi, omega, m, rep=100){

  p <- ncol(phi)
  truth <- target.sigma(phi, omega)
  start <- matrix(0, nrow = m, ncol = p)  #only depends on C

  for(i in floor(m/2):1){
    start[i,] <- 2*i*sqrt(diag(truth))
    start[m-i+1,] <- -i*sqrt(diag(truth))
  }
  master.chain.rep <- array(0,dim = c(max, p, m, rep))

  for (r in 1:rep){
    for (k in 1:m){
      master.chain.rep[,,k,r] <- markov.chain(phi, omega, max, start[k,])
    }
  }

  conv.pts <- seq(min, max,100)
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
  save(asv.samp,rsv.samp, ess.asv.samp, ess.rsv.samp, file = paste(paste("Out/bartlett/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

}

#########################################
###Tukey-Hanning
#########################################

create.output2 <- function(phi, omega, m, check.pts, freq, c.prob){
  
  p <- ncol(phi)
  truth <- target.sigma(phi, omega)
  start <- matrix(0, nrow = m, ncol = p)  #only depends on C
  
  for(i in floor(m/2):1){
    start[i,] <- 2*i*sqrt(diag(truth))
    start[m-i+1,] <- -2*i*sqrt(diag(truth))
  }
  
  critical <- qchisq(c.prob, df=2)
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    
    asv.samp <- array(0, dim = c(p,p,freq))
    rsv.samp <- array(0, dim = c(p,p,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    for (j in 1:freq){
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      chain <- array(0,dim = c(nsim,p,m))
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        chain[,,k] <- markov.chain(phi, omega, nsim, start[k,])
        #b[k] <- batchSize(chain[,,k], method = "bartlett")
        model <- lm(chain[,,k] ~ 1)
        b[k] <- bwAndrews(model, kernel = "Tukey-Hanning", approx = "AR(1)", prewhite = 0)
      }
      
      b.avg <- floor(mean(b))
      global.mean <- apply(chain, 2, mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg, method = "tukey")
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg, method = "tukey")
      }
      
      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
      
      chi.sq.asv <- t2.stat(global.mean, rep(0,p), asv.samp[,,j], nsim*m)
      chi.sq.rsv <- t2.stat(global.mean, rep(0,p), rsv.samp[,,j], nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp,asv.coverage,rsv.coverage, file = paste(paste("Out/tukey/out", m, nsim, sep = "_"),".Rdata", sep = ""))
  }
  
}

##################################################
#function for storing convergence plots data
##################################################
convergence2 <- function(min, max, phi, omega, m, rep=100){
  
  p <- ncol(phi)
  truth <- target.sigma(phi, omega)
  start <- matrix(0, nrow = m, ncol = p)  #only depends on C
  
  for(i in floor(m/2):1){
    start[i,] <- 2*i*sqrt(diag(truth))
    start[m-i+1,] <- -i*sqrt(diag(truth))
  }
  master.chain.rep <- array(0,dim = c(max, p, m, rep))
  
  for (r in 1:rep){
    for (k in 1:m){
      master.chain.rep[,,k,r] <- markov.chain(phi, omega, max, start[k,])
    }
  }
  
  conv.pts <- seq(min, max,100)
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
        #b[k] <- batchSize(chain[,,k], method = "bartlett")
        model <- lm(chain[,,k] ~ 1)
        b[k] <- bwAndrews(model, kernel = "Tukey-Hanning", approx = "AR(1)", prewhite = 0)
        smpl.cov[,,k] <- cov(chain[,,k])
      }
      
      b.avg <- floor(mean(b))
      global.mean <- apply(chain,2,mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg, method = "tukey")
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg, method = "tukey")
        
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
  save(asv.samp,rsv.samp, ess.asv.samp, ess.rsv.samp, file = paste(paste("Out/tukey/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))
  
}


######################################################################

p <- 5
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .5^i
    omega[j+i, j] <- .5^i
  }
}

phi <- diag(c(.99, .1, .1, .1, .1))

m <- 2

#sims for plotting densities and calculating coverage

check.pts <- c(3e2, 5e2, 1e3, 5e3, 1e4)
freq <- 1e3  #100 for now, will change later
rep <- 10
c.prob <- .95
min <- 5e2
max <- 1e4
conv.pts <- seq(min, max, 100)

print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output1(phi, omega, m, check.pts, freq, c.prob)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")
convergence1(min, max, phi, omega, m, rep)

print("Carrying out 1000 repititions for each value of nsim in check.pts, method = tukey")
create.output2(phi, omega, m, check.pts, freq, c.prob)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5, method = tukey")
convergence2(min, max, phi, omega, m, rep)

