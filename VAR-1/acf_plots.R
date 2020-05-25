set.seed(4)
source("functions.R")
library(rep.acf.ccf)


m <- 2
nsim <- 5e2
p <- 2
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .5^i
    omega[j+i, j] <- .5^i
  }
}

p <- 2
phi <- diag(c(.9999, .001))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

# phi <- phi/(max(eigen(phi)$values) + .1)

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

for(i in 1:floor(m/2)){
  start[i,] <- 1*i*sqrt(diag(sigma))
  start[m-i+1,] <- -1*i*sqrt(diag(sigma))
}

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m
plot(mc.chain.list[[1]], xlim = range(c(mc.chain.list[[1]][,1],mc.chain.list[[2]][,1]) ), ylim = range(c(mc.chain.list[[1]][,2],mc.chain.list[[2]][,2]) ), asp = 1)
points(mc.chain.list[[2]], col = "red")
points(mc.chain.list[[3]], col = "blue")

plot.ts(mc.chain.list[[1]])
plot.ts(mc.chain.list[[2]])
plot.ts(mc.chain.list[[3]])

acf.list <- combined_acf(mc.chain.list, chain = 1, component = seq(1,p), lag.max = 150, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "correlation")

# pdf(file = paste("Out/acf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim))
par(mfrow = c(3,2))

plot(acf.list[[1]][[1]], main = paste("Locally centered ACF plot for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[1]][[2]], main = paste("Globally centered ACF plot for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[2]][[1]], main = paste("Locally centered ACF plot for component - 5"), xlim = range(0, lag.max),ylim =  c(0,1))#c(min(min(acf.list[[2]][[1]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[1]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")
plot(acf.list[[2]][[2]], main = paste("Globally centered ACF plot for component - 5"), xlim = range(0, lag.max), ylim =  c(0,1)) #c(min(min(acf.list[[2]][[2]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[2]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")

plot(ccf.list[[1]], main = "Locally centered CCF plot between 1 and 5", xlim = c(-lag.max, lag.max), ylim = range(ccf.list[[1]]$acf))
lines(seq(-lag.max, lag.max), true.acf[1,2,]/true.acf[1,2,lag.max+1], col = "red")
plot(ccf.list[[2]], main = "Globally centered CCF plot between 1 and 5", ylim = range(ccf.list[[2]]$acf))
lines(seq(-lag.max, lag.max), true.acf[1,2,]/true.acf[1,2,lag.max+1], col = "red")
# dev.off()
