set.seed(1)
source("functions.R")
library(rep.acf.ccf)

m <- 5
p <- 2
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(c(.9999, .001))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)
#phi <- phi/(max(eigen(phi)$values) + .1)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)

lag.max <- 150
true.acf <- array(0, dim = c(p, p, 2*lag.max + 1))
true.acf[,,lag.max+1] <- target
for (i in 1:lag.max){
  true.acf[,,lag.max + 1 + i] <- phi %*% true.acf[,,lag.max + 1 + i - 1]
  true.acf[,,lag.max + 1 - i] <- true.acf[,,lag.max + 1 - i + 1] %*% t(phi)
}

start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- 0.5*i*sqrt(diag(target))
  start[m-i+1,] <- -0.5*i*sqrt(diag(target))
}

####################################
nsim <- 500
########################################
mc.chain.list <- list()
global.mean <- rep(0,p)
lag.max <- 150
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
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

pdf(file = "acf_comparisons_5e2.pdf", width = 6, height = 4)
plot(acf.list[[1]][[1]], main = expression(paste("ACF comparison for chain 1")), type = "l", lty = 1, xlim = c(0,lag.max), ylim = c(-0.5,1))
lines(seq(0,lag.max), acf.list[[1]][[2]]$acf, type = "l", col = "blue", lty = 2, xlim = c(0,lag.max), ylim = c(0,1))
lines(seq(0, lag.max), new_corr, col = "dark green", lty = 4)
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red", lty = 5)
legend("topright", legend=c("ACF", "R-ACF", "Split", "Truth"),col=c("black", "blue", "dark green", "red"), lty=c(1,2,4,5), cex=.75)
dev.off()


####################################
nsim <- 10000
########################################

mc.chain.list <- list()
global.mean <- rep(0,p)
lag.max <- 150
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
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

pdf(file = "acf_comparisons_1e4.pdf", width = 6, height = 4)
plot(acf.list[[1]][[1]], main = expression(paste("ACF comparison for chain 1")), type = "l", lty = 1, xlim = c(0,lag.max), ylim = c(-0.5,1))
lines(seq(0,lag.max), acf.list[[1]][[2]]$acf, type = "l", col = "blue", lty = 2, xlim = c(0,lag.max), ylim = c(0,1))
lines(seq(0, lag.max), new_corr, col = "dark green", lty = 4)
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red", lty = 5)
legend("topright", legend=c("ACF", "R-ACF", "Split", "Truth"),col=c("black", "blue", "dark green", "red"), lty=c(1,2,4,5), cex=.75)
dev.off()

