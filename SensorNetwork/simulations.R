set.seed(10)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")
######## Data

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

############################################################
##creates freq=1000 replications of ASV and RSV for each value of nsim from check.pts
############################################################

create.output <- function(start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, check.pts, freq,truth){

  p <- 8

  for (i in 1:length(check.pts)){

    nsim <- check.pts[i]
    critical <- ((nsim*m - 1)*p/(nsim*m  - p))*qf(.95, df1 = p, df2 = (nsim*m-p))
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

      for(k in 1:m){
        temp <- MHwG.RAM(start[k,], aux[k,], jump.scale = j.scale, Ob, Os, Xb, Xs, Yb, Ys,
                               n.sample = nsim+1, n.burn = 0)
        chain[,,k] <- temp$x
        print(colMeans(temp$accept))
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

      chi.sq.asv <- t2.stat(global.mean,truth,asv.samp[,,j],nsim*m)
      chi.sq.rsv <- t2.stat(global.mean,truth,rsv.samp[,,j],nsim*m)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.coverage, rsv.coverage, asv.samp,rsv.samp, file = paste(paste("Out/out",nsim, sep = "_"),".Rdata", sep = ""))
  }

}

##################################################
#function for storing convergence plots data
##################################################
convergence <- function(min, max, start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, rep=100, step=500){

  p <- 8

  #master.chain.rep <- array(0, dim = c(max, p, m, rep))

  #for (i in 1:rep){
   # print(paste("Sampling for rep =", i))
    #for (j in 1:m){
     # print(paste("Sampling chain = ", j))
    # temp <- MHwG.RAM(start[j,], aux[j,], jump.scale = j.scale,Ob, Os, Xb, Xs, Yb, Ys, n.sample = max+1, n.burn = 0)
    #  master.chain.rep[,,j,i] <- temp$x
    #  print(colMeans(temp$accept))
    #}
  #}

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
      
      print(paste("Sampling for nsim = ", nsim, "; rep = ", r))
      
      chain <-array(0, dim = c(nsim,p,m))
      #chain <- master.chain.rep[1:nsim,,,r]
      chain[,,1] <- MHwG.RAM(start[1,], aux[1,], jump.scale = j.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample = nsim+1, n.burn = 0)$x
      chain[,,2] <- MHwG.RAM(start[2,], aux[2,], jump.scale = j.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample = nsim+1, n.burn = 0)$x
      
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
  save(asv.samp,rsv.samp, ess.asv.samp, ess.rsv.samp, file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

}


######################################################################
m <- 2
start1 <- runif(n=8, min=-0.3, max=0)
aux1 <- runif(n=8, min=-0.3, max=0)
start2 <- runif(n=8, min= 0.7, max=1)
aux2 <- runif(n=8, min= 0.7, max=1)
start <- rbind(start1, start2)
aux <- rbind(aux1, aux2)
j.scale <- rep(1,.08, 4)
truth <- c(0.5748, 0.9069, 0.0991, 0.3651, 0.2578, 0.1350, 0.8546, 0.0392)
#sims for plotting densities and calculating coverage

check.pts <- c(1e3, 5e3, 1e4, 5e4, 1e5)
freq <- 10
rep <- 2
c.prob <- .95
min <- 1e3
max <- 2e5
step <- 5e2
conv.pts <- seq(min, max, step)

print("Carrying out 1000 repititions for each value of nsim in check.pts")
create.output(start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, check.pts, freq, truth)

print("Carrying out simulations for convergence plots of ASV and RSV in the range(1e3, 1e5")

convergence(min, max, start, aux, j.scale, Ob, Os, Xb, Xs, Yb, Ys, m, rep, step=500)


