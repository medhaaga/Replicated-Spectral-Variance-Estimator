################################################################
##This code is responsible for simulating and storing R objects
##that are later used for plotting figures and making tables.
################################################################

set.seed(10)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")
sourceCpp("MCMC.cpp")


##################################################
#function for storing convergence plots data
##################################################

convergence <- function(min, max, step, phi, omega, target, m, rep=100, set){

  p <- ncol(phi)
  start <- rmvnorm(m, mean = rep(0, p), sigma = target)

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
      master.chain[,,k] <- markov_chain(phi, omega, max, start[k,])
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
      ess.asv.samp[j] <- exp((1/p)*(sum(log(eigen(lambda.rep)$values)) - sum(log(eigen(asv.samp[,,j])$values))))
      ess.rsv.samp[j] <- exp((1/p)*(sum(log(eigen(lambda.rep)$values)) - sum(log(eigen(rsv.samp[,,j])$values))))

    }

    asv[[r]] <- asv.samp
    rsv[[r]] <- rsv.samp
    ess.asv[[r]] <- ess.asv.samp
    ess.rsv[[r]] <- ess.rsv.samp

    save(asv, rsv, ess.asv, ess.rsv,
         file = paste("Out/conv_data_min", min, "_max", max, "_set", set, ".Rdata", sep = ""))
  }

}

######################################################################

m <- 5 #### number of chains
p <- 100  #### dimension


##################################################
####### Setting 1: Good mixing ###################
##################################################

#### target distribution specifications######

set <- 1
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(0.1 + seq(0, (p-1))*((0.8 - 0.1)/(p-1)))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
save(target, truth, file = "Out/var-set1_truth.Rdata")

rep <- 10
min <- 5e2
max <- 5e4
step <- 500
conv.pts <- seq(min, max, step)

start.time <- Sys.time()
print("Carrying out simulations for convergence plots of ASV and RSV in the range(5e2, 5e4) for good mixing")
convergence(min, max, step, phi, omega, target, m, rep, set)
end.time <- Sys.time()
print(end.time - start.time)


##################################################
####### Setting 2: Bad mixing ###################
##################################################

#### target distribution specifications######

set <- 2
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(0.99 + seq(0, (p-1))*((0.999 - 0.99)/(p-1)))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
save(target, truth, file = "Out/var-set2_truth.Rdata")

rep <- 10
min <- 5e2
max <- 5e4
step <- 500
conv.pts <- seq(min, max, step)

start.time <- Sys.time()
print("Carrying out simulations for convergence plots of ASV and RSV in the range(5e2, 5e4) for bad mixing")
convergence(min, max, step, phi, omega, target, m, rep, set)
end.time <- Sys.time()
print(end.time - start.time)

