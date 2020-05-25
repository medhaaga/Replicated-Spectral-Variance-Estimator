set.seed(40)
source("functions.R")
library(rep.acf.ccf)


m <- 2
nsim <- 1000
p <- 5
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .5^i
    omega[j+i, j] <- .5^i
  }
}


phi <- diag(c(.99, .1, .1, .1, .1))
sigma <- target.sigma(phi, omega)
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

for(i in floor(m/2):1){
  start[i,] <- 2*i*sqrt(diag(sigma))
  start[m-i+1,] <- -2*i*sqrt(diag(sigma))
}

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, chain = 1, component = seq(1,p), lag.max = 150, type = "covariance")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,5), lag.max = 150, type = "covariance")

pdf(file = paste("Out/acf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim))
par(mfrow = c(3,2))

plot(acf.list[[1]][[1]], main = paste("Locally centered ACF plot for component - 1"), xlim = range(0, lag.max), ylim = c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,], col = "red")
plot(acf.list[[1]][[2]], main = paste("Globally centered ACF plot for component - 1"), xlim = range(0, lag.max), ylim = c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,], col = "red")
plot(acf.list[[5]][[1]], main = paste("Locally centered ACF plot for component - 5"), xlim = range(0, lag.max), ylim = c(min(min(acf.list[[5]][[1]]$acf), min(true.acf[5,5,])), max(max(acf.list[[5]][[1]]$acf), max(true.acf[5,5,]))))
lines(seq(-lag.max, lag.max), true.acf[5,5,], col = "red")
plot(acf.list[[5]][[2]], main = paste("Globally centered ACF plot for component - 5"), xlim = range(0, lag.max), ylim = c(min(min(acf.list[[5]][[2]]$acf), min(true.acf[5,5,])), max(max(acf.list[[5]][[2]]$acf), max(true.acf[5,5,]))))
lines(seq(-lag.max, lag.max), true.acf[5,5,], col = "red")

plot(ccf.list[[1]], main = "Locally centered CCF plot between 1 and 5", xlim = c(-lag.max, lag.max), ylim = range(ccf.list[[1]]$acf))
lines(seq(-lag.max, lag.max), true.acf[1,5,], col = "red")
plot(ccf.list[[2]], main = "Globally centered CCF plot between 1 and 5", ylim = range(ccf.list[[2]]$acf))
lines(seq(-lag.max, lag.max), true.acf[1,5,], col = "red")
dev.off()
