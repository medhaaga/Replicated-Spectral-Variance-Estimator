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


##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, rep=100, step=500){

  p <- 8
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
      temp <- MHwG.RAM(start[k,], aux[k,], jump.scale = j.scale,Ob, Os, Xb, Xs, Yb, Ys, n.sample = max+1, n.burn = 0)
      master.chain[,,k] <- temp$x
      print(colMeans(temp$accept)) 
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
    file = paste("Out/conv_data_m", m, "_min", min, "_max", max, ".Rdata", sep = "")) 
  }
  
}

#############################################
################### Data#####################
#############################################

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors.
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)

# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546,
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0,
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)

####################################################
############ Model Parameters#######################
###################################################

m <- 5 #### number of chains
p <- 8 #### dimensions

####### Starting Values####################

start1 <- c(-0.1, 0.5, -0.1, -0.2, 0.1, 0.1, -0.5, -0.5)
aux1 <- runif(n=8, min=min(start1), max=max(start1))
start2 <- c(0.0, 0.6, 0.1, 0.1, 0.2, 0.2, 1.0, 0.0)
aux2 <- runif(n=8, min=min(start2), max=max(start2))
start3 <- c(0.2, 0.7, 0.5, 0.4, 0.5, 0.3, 0.5, 0.5)
aux3 <- runif(n=8, min=min(start3), max=max(start3))
start4 <- c(0.4, 0.8, 0.8, 0.6, 0.7, 0.4, 1.0, 1.0)
aux4 <- runif(n=8, min=min(start4), max=max(start4))
start5 <- c(0.7, 1.0, 1.2, 0.9, 0.9, 0.5, 1.5, 1.5)
aux5 <- runif(n=8, min=min(start5), max=max(start5))

start <- rbind(start1, start2, start3, start4, start5)
aux <- rbind(aux1, aux2, aux3, aux4, aux5)
j.scale <- rep(0.5, 4)

#sims for plotting running plots

rep <- 100
min <- 500
max <- 2e5
step <- 500
conv.pts <- seq(min, max, step)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")
convergence(min, max, start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, rep, step=500)


