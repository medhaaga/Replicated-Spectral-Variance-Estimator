set.seed(1)
source("functions.R")
library(rep.acf.ccf)


m <- 5
A <- 1
B <- 3
C <- 8
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C


for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

####################################
nsim <- 500
########################################
mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")

pdf(file = paste(paste("Out/", A, B, C, "/acf, n=", sep = "_"), nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim), height = 3)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = paste("ACF for component - 1"))
plot(acf.list[[1]][[2]], main = paste("R-ACF for component - 1"))

dev.off()

####################################
nsim <- 1000
########################################
mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")

pdf(file = paste(paste("Out/", A, B, C, "/acf, n=", sep = "_"), nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim), height = 3)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = paste("ACF for component - 1"))
plot(acf.list[[1]][[2]], main = paste("R-ACF for component - 1"))

dev.off()

####################################
nsim <- 10000
########################################
mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "covariance")

pdf(file = paste(paste("Out/", A, B, C, "/acf, n=", sep = "_"), nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim), height = 3)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = paste("ACF for component - 1"))
plot(acf.list[[1]][[2]], main = paste("R-ACF for component - 1"))

dev.off()


