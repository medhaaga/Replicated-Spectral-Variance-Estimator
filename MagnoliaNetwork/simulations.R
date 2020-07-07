set.seed(1)

library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")
library(ergm)

data(faux.magnolia.high)
magnolia = faux.magnolia.high

#Ensure well connected component only 
notwellconnected = which(component.largest(magnolia, connected="weak")==FALSE)
delete.vertices(magnolia, notwellconnected)
dp = dataPrepHighSchool(magnolia)

############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output <- function(start, m, check.pts, freq, truth){
  
  p <- 5
  
  for (i in 1:length(check.pts)){
    
    nsim <- check.pts[i]
    critical <- ((nsim*m - 1)*p/(nsim*m  - p))*qf(.95, df1 = p, df2 = (nsim*m-p))
    asv.samp <- array(0, dim = c(p,p,freq))
    rsv.samp <- array(0, dim = c(p,p,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    
    for (j in 1:freq){
      
      print(j)
      if(j %% (freq/10) == 0) print(paste("Percentage completion: ", round(j/freq*100, 2), "for nsim = ", nsim))
      
      chain <- array(0,dim = c(nsim,p,m))
      sve <- array(0, dim = c(p,p,m))
      rsve <- array(0, dim = c(p,p,m))
      b <- rep(0,m)
      
      temp <- magnolia.mhrw.func(nrep = nsim, chains = m, start = start)
      for (k in 1:m){
        chain[,,k] <- as.matrix(MHRWestimation(temp[[k]])[,c(1,7,5,8,9)], nrow = nsim, ncol = p)
        b[k] <- batchSize(chain[,,k], method = "bartlett")
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
convergence <- function(min, max, start, m, rep=100, step=500){
  
  p <- 5  #dimension of sample space
  
  master.chain.rep <- array(0, dim = c(max, p, m, rep))
  
  for (i in 1:rep){
    print(paste("Sampling for rep =", i))
    for (j in 1:m){
      print(paste("Sampling chain = ", j))
      temp <- MHRWestimation(magnolia.mhrw.func(nrep = max, chains = 1, start = start[j])[[1]])
      master.chain.rep[,,j,i] = as.matrix(temp[,c(1,7,5,8,9)])
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
    
      chain <- master.chain.rep[1:nsim,,,r]
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
      ess.asv.rep[r] <- (det(lambda.rep[,,r])/det(sve.rep[,,r]))^(1/p) #(det(lambda.rep[,,r])/det(sve.rep[,,r]))^(1/p)
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
m <- 2
start <- c(279, 360)

#sims for plotting densities and calculating coverage

check.pts <- c(1e2, 5e2, 1e3, 5e3, 1e4)
freq <- 20
rep <- 10
c.prob <- .95
min <- 1e2
max <- 1e4
step <- 1e2
conv.pts <- seq(min, max, step)
truth = c(meandegree = mean(dp$degree.connected$degree),
             meancluster = mean(ifelse(dp$degree.connected$degree>1, dp$degree.connected$triangles/choose(dp$degree.connected$degree,2), 0)),
             meangrade = mean(dp$degree.connected$grade), 
             propfemale = mean(ifelse(dp$degree.connected$sex=="F", 1, 0)), 
             propwhite = mean(ifelse(dp$degree.connected$race=="White", 1, 0)))

print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output(start, m, check.pts, freq, truth)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")

convergence(min, max, start, m, rep, step)

