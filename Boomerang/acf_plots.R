set.seed(1)
source("functions.R")
library(rep.acf.ccf)


m <- 2
A <- 1
B <- 3
C <- 8
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C


for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

lag.max <- 100
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

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf,n=", sep = "_"), nsim, ".pdf", sep = ""), width = 8, height = 4)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = expression(paste("Old ACF for chain 1")))
plot(acf.list[[1]][[2]], main = expression(paste("New ACF for chain 1")))

dev.off()

pdf(file = paste(paste("Out/", A, B, C, "/scatter_n", sep = "_"), nsim, ".pdf", sep = ""), width = 5, height = 5)
par(mfrow = c(1,1))
plot(mc.chain.list[[1]], xlim = c(0,11), ylim = c(0,11), col = 1, pch = 1)
points(mc.chain.list[[2]], col = "red", pch = 2)
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

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf,n=", sep = "_"), nsim, ".pdf", sep = ""), width = 8, height = 4)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = expression(paste("Old ACF for chain 1")))
plot(acf.list[[1]][[2]], main = expression(paste("New ACF for chain 1")))

dev.off()


pdf(file = paste(paste("Out/", A, B, C, "/scatter_n", sep = "_"), nsim, ".pdf", sep = ""), width = 5, height = 5)
par(mfrow = c(1,1))
plot(mc.chain.list[[1]], xlim = c(0,11), ylim = c(0,11), col = 1, pch = 1)
points(mc.chain.list[[2]], col = "red", pch = 2)
dev.off()
####################################
nsim <- 50000
########################################
mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf,n=", sep = "_"), nsim, ".pdf", sep = ""), width = 8, height = 4)
#par(mfrow = c(3,2))
#for (i in 1:2){
#  plot(acf.list[[i]][[1]], main = paste("ACF for component -", i ))
#  plot(acf.list[[i]][[2]], main = paste("R-ACF for component -", i ))
#}
#plot(ccf.list[[1]], main = "CCF for component 1-2")
#plot(ccf.list[[2]], main = "R-CCF for component 1-2")
par(mfrow = c(1,2))
plot(acf.list[[1]][[1]], main = expression(paste("Old ACF for chain 1")))
plot(acf.list[[1]][[2]], main = expression(paste("New ACF for chain 1")))

dev.off()

pdf(file = paste(paste("Out/", A, B, C, "/scatter_n", sep = "_"), nsim, ".pdf", sep = ""), width = 5, height = 5)
par(mfrow = c(1,1))
plot(mc.chain.list[[1]], xlim = c(0,11), ylim = c(0,11), col = 1, pch = 1)
points(mc.chain.list[[2]], col = "red", pch = 2)
dev.off()


