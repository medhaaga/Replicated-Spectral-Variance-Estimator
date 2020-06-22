set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(spBayes)
source("functions.R")
sourceCpp("lag.cpp")
load("spatio_temp_truth_2e6by1000")
data("NETemp.dat")
ne.temp <- NETemp.dat



############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output <- function(N.t, max.d, m, check.pts, freq, c.prob, truth){
  
  p <- 185
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    critical <- ((nsim*m-1)*p/(nsim*m  - p))*qf(.95, df1 = p, df2 = (nsim*m-p))
    asv.samp <- array(0, dim = c(p,p,freq))
    rsv.samp <- array(0, dim = c(p,p,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    starting1 <- list("beta"=rep(30,N.t*2), "phi"=rep(3/(0.1*max.d), N.t),
                      "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                      "sigma.eta"=diag(rep(0.01, 2)))
    
    starting2 <- list("beta"=rep(-30,N.t*2), "phi"=rep(3/(0.85*max.d), N.t),
                      "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                      "sigma.eta"=diag(rep(0.01, 2)))
    
    for (j in 1:freq){
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      chain <- array(0,dim = c(nsim,p,m))
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      chain[,,1] <- spatio_temp(N = nsim, verbose = FALSE, starting1)
      chain[,,2] <- spatio_temp(N = nsim, verbose = FALSE, starting2)
      b[1] <- batchSize(chain[,,1], method = "bartlett")
      b[2] <- batchSize(chain[,,2], method = "bartlett")
      
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
      
      chi.sq.asv <- t2.stat(global.mean,truth,asv.samp[,,j],nsim*m)
      chi.sq.rsv <- t2.stat(global.mean,truth,rsv.samp[,,j],nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp, asv.coverage, rsv.coverage, file = paste(paste("Out/out",nsim, sep = "_"),".Rdata", sep = ""))
  }
  
}

##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, N.t, max.d, m, rep=100, step=500){
  
  p <- 185
  
  starting1 <- list("beta"=rep(30,N.t*2), "phi"=rep(3/(0.1*max.d), N.t),
                    "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                    "sigma.eta"=diag(rep(0.01, 2)))
  
  starting2 <- list("beta"=rep(-30,N.t*2), "phi"=rep(3/(0.8*max.d), N.t),
                    "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                    "sigma.eta"=diag(rep(0.01, 2)))
 
  
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
      chain <-array(0, dim = c(nsim,p,m))
      chain[,,1] <- spatio_temp(N = nsim, verbose = FALSE, starting1)
      chain[,,2] <- spatio_temp(N = nsim, verbose = FALSE, starting2)
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


######################################################################
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]##subset first 2 years (Jan 2000 - Dec. 2002)

coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))
N.t <- 12
m <- 2

#sims for plotting densities and calculating coverage

check.pts <- c(5e2, 1e3, 5e3, 1e4)
freq <- 1e2  
rep <- 10
c.prob <- .95
min <- 5e2
max <- 5e4
step <- 100
conv.pts <- seq(min, max, step)
truth <- colMeans(spatio_truth)


print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output(N.t, max.d, m, check.pts, freq, c.prob, truth)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")

convergence(min, max, N.t, max.d, m, rep, step)

