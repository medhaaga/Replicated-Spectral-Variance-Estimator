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
sample_var <- rep(0, m)
means <- rep(0, m)
covariances <- matrix(0, nrow = lag.max+1, ncol = m)

for (j in 1:m){
  sample_var[j] <- var(mc.chain.list[[j]][,component])
  means[j] <- mean(mc.chain.list[[j]][,component])
  covariances[,j] <- as.vector(acf(mc.chain.list[[j]][,component], plot = FALSE, lag.max = lag.max, type = "covariance")$acf)
}
W <- mean(sample_var)
B <- nsim * var(means)
new_var <- ((nsim - 1)/nsim)*W + (1/nsim)*B
avg.covariances <- as.matrix(rowMeans(covariances))
new_corr <- 1 - (W/new_var) + (avg.covariances/new_var)

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = "Boomerang.pdf", width = 8, height = 4)
plot(acf.list[[1]][[1]], main = expression(paste("ACF comparison for chain 1")), type = "l", xlim = c(0,lag.max), ylim = c(-0.5,1))
lines(seq(0,lag.max), acf.list[[1]][[2]]$acf, type = "l", col = "blue", xlim = c(0,lag.max), ylim = c(0,1))
lines(seq(0, lag.max), new_corr, col = "green")
legend("topright", legend=c("ACF", "R-ACF", "Split"),col=c("black", "blue", "green"), lty=1, cex=.75)
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

for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m
sample_var <- rep(0, m)
means <- rep(0, m)
covariances <- matrix(0, nrow = lag.max+1, ncol = m)

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf,n=", sep = "_"), nsim, ".pdf", sep = ""), width = 8, height = 4)
plot(acf.list[[1]][[1]], main = expression(paste("ACF comparison for chain 1")), type = "l", xlim = c(0,lag.max), ylim = c(-0.5,1))
lines(seq(0,lag.max), acf.list[[1]][[2]]$acf, type = "l", col = "blue", xlim = c(0,lag.max), ylim = c(0,1))
lines(seq(0, lag.max), new_corr, col = "green")
legend("topright", legend=c("ACF", "R-ACF", "Split"),col=c("black", "blue", "green"), lty=1, cex=.75)
dev.off()



####################################
nsim <- 50000
####################################

mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m
sample_var <- rep(0, m)
means <- rep(0, m)
covariances <- matrix(0, nrow = lag.max+1, ncol = m)

acf.list <- combined_acf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = lag.max, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf,n=", sep = "_"), nsim, ".pdf", sep = ""), width = 8, height = 4)
plot(acf.list[[1]][[1]], main = expression(paste("ACF comparison for chain 1")), type = "l", xlim = c(0,lag.max), ylim = c(-0.5,1))
lines(seq(0,lag.max), acf.list[[1]][[2]]$acf, type = "l", col = "blue", xlim = c(0,lag.max), ylim = c(0,1))
lines(seq(0, lag.max), new_corr, col = "green")
legend("topright", legend=c("ACF", "R-ACF", "Split"),col=c("black", "blue", "green"), lty=1, cex=.75)
dev.off()



