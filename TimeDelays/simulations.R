set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(timedelay)
source("functions.R")
sourceCpp("lag.cpp")
data <- read.csv("q0957usno.csv", header = TRUE)
lcA <- data[, 1 : 3]
lcB <- data[, c(1, 4, 5)]

################################################3

############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output <- function(lcA, lcB, micro, ram, delta.start, delta.jump, m, check.pts, freq){
  
  p <- 2
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    #critical <- ((nsim*m - 1)*p/(nsim*m  - p))*qf(.95, df1 = p, df2 = (nsim*m-p))
    asv.samp <- array(0, dim = c(p,p,freq))
    rsv.samp <- array(0, dim = c(p,p,freq))
    #asv.coverage <- rep(0,freq)
    #rsv.coverage <- rep(0,freq)
    
    
    for (j in 1:freq){
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      chain <- array(0,dim = c(nsim,p,m))
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for(k in 1:m){
        chain[,,k] <- markov.chain(lcA, lcB, nsim+1, delta.start[k], delta.jump, micro, ram, burn = 0)$samples
        b[k] <- batchSize(chain[,,k], method = "bartlett")
        }
      
      b.avg <- floor(mean(b))
      global.mean <- apply(chain,2,mean)
      
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
      }
      
      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
      
      #chi.sq.asv <- t2.stat(global.mean,truth,asv.samp[,,j],nsim*m)
      #chi.sq.rsv <- t2.stat(global.mean,truth,rsv.samp[,,j],nsim*m)
      #if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      #if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp, file = paste(paste("Out2/out",nsim, sep = "_"),".Rdata", sep = ""))
  }
  
}

##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, lcA, lcB, micro, delta.start, delta.jump, ram, m, rep=100, step=500){
  
  p <- 2
  
  master.chain.rep <- array(0, dim = c(max, p, m, rep))
  for (i in 1:rep){
    print(paste("Sampling for rep =", i))
    for (j in 1:m){
      print(paste("Sampling chain = ", j))
      master.chain.rep[,,j,i] <- markov.chain(lcA, lcB, max+1, delta.start[j], delta.jump = 100, micro, ram, burn = 0)$samples
    }
  }
  
  conv.pts <- seq(min, max, step)
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
      #chain <-array(0, dim = c(nsim,p,m))
      chain <- master.chain.rep[1:nsim,,,r]
      #chain[,,1] <- spatio_temp(N = nsim, verbose = FALSE, starting1)
      #chain[,,2] <- spatio_temp(N = nsim, verbose = FALSE, starting2)
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      smpl.cov <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      for (k in 1:m){
        b[k] <- batchSize(chain[,,k], method = "bartlett")
        smpl.cov[,,k] <- var(chain[,,k])
      }
      
      b.avg <- floor(mean(b))
      
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
      ess.asv.rep[r] <- (lambda.rep[1,1,r]/sve.rep[1,1,r]) #(det(lambda.rep[,,r])/det(sve.rep[,,r]))^(1/p)
      ess.rsv.rep[r] <- (lambda.rep[1,1,r]/rsve.rep[1,1,r])
    }
    
    asv.samp[,,j] <- apply(sve.rep, c(1,2), mean)
    rsv.samp[,,j] <- apply(rsve.rep, c(1,2), mean)
    ess.asv.samp[[j]] <- ess.asv.rep
    ess.rsv.samp[[j]] <- ess.rsv.rep
    print(paste("Percentage completion: ", round(100*j/l, 2)))
  }
  save(asv.samp,rsv.samp, ess.asv.samp, ess.rsv.samp, file = paste(paste("Out2/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))
  
}


######################################################################
m <- 4
delta.start <- c(300, 500, 700, 900)
delta.jump <- 100
micro <- 0
ram <- TRUE

#sims for plotting densities and calculating coverage

check.pts <- c(5e2, 1e3, 5e3, 1e4, 5e4)
freq <- 100 
rep <- 10
c.prob <- .95
min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)
#truth <- c(50, 2, (colMeans(chain[,,1])[3:158] + colMeans(chain[,,2])[3:158])/2, 0, 0.03, 100)

print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output(lcA, lcB, micro, ram, delta.start, delta.jump, m, check.pts, freq)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")

convergence(min, max, lcA, lcB, micro, delta.start, delta.jump, ram, m, rep, step=500)

