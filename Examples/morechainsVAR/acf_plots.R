set.seed(10)
source("functions.R")
library(multichainACF)
library("mvtnorm")

######################################
##### Model parameters ###############
######################################

m <- 101
p <- 2
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(c(.999, .001))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
lag.max <- 40

true.acf <- array(0, dim = c(p, p, 2*lag.max + 1))
true.acf[,,lag.max+1] <- target
for (i in 1:lag.max){
  true.acf[,,lag.max + 1 + i] <- phi %*% true.acf[,,lag.max + 1 + i - 1]
  true.acf[,,lag.max + 1 - i] <- true.acf[,,lag.max + 1 - i + 1] %*% t(phi)
}

################ Only for creating Markov chains. Don't run. ##################
start <- matrix(0, nrow = m, ncol = p)  #only depends on C
for(i in 1:m)
  start[i,] <- rmvnorm(n=1, mean = rep(0,p), sigma = target)
# for(i in 1:floor(m/2)){
#   start[(2*i),] <- i*.5*sqrt(diag(target))
#   start[(2*i +1),] <- -i*.5*sqrt(diag(target))
# }

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  print(i)
  print(paste("Starting at: ", start[i,]))
  chain <- markov.chain(phi, omega, 1e4, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
  print(paste("Local mean: ", colMeans(chain)))
}
global.mean <- global.mean/m


save(mc.chain.list, true.acf, file = "Out/var-100_chains.Rdata")
#################################################################################
