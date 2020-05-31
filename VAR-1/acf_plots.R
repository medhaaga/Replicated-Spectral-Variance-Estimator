set.seed(10)
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

##################################
nsim <- 5e2
##################################

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

pdf(file = paste("Out/scatter_plot_", nsim, ".pdf", sep = ""), title = paste("Scatter plot, nsim = ", nsim))
plot(mc.chain.list[[1]], xlim = c(-(m/2)*sqrt(diag(target)[1]), (m/2)*sqrt(diag(target)[1])), ylim = c(-(m/2)*sqrt(diag(target)[2]), (m/2)*sqrt(diag(target)[2])), xlab = "X component", ylab = "Y component")
points(mc.chain.list[[2]], col = 2)
points(mc.chain.list[[3]], col = 3)
points(mc.chain.list[[4]], col = 4)
points(mc.chain.list[[5]], col = 5)
dev.off()

acf.list <- combined_acf(mc.chain.list, chain = 2, component = seq(1,p), lag.max = 150, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 2, component = c(1,2), lag.max = 150, type = "correlation")

pdf(file = paste("Out/acf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim), height = 4)
par(mfrow = c(2,2))

plot(acf.list[[1]][[1]], main = paste("ACF for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[1]][[2]], main = paste("R-ACF for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[2]][[1]], main = paste("ACF for component - 2"), xlim = range(0, lag.max),ylim =  c(0,1))#c(min(min(acf.list[[2]][[1]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[1]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")
plot(acf.list[[2]][[2]], main = paste("R-ACF for component - 2"), xlim = range(0, lag.max), ylim =  c(0,1)) #c(min(min(acf.list[[2]][[2]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[2]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")
dev.off()


pdf(file = paste("Out/ccf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered CCF plot, nsim = ", nsim), height = 3)
par(mfrow = c(1,2))
plot(ccf.list[[1]], main = "CCF between 1 and 2", xlim = c(-lag.max, lag.max), ylim = range(ccf.list[[1]]$acf))
lines(seq(-lag.max, lag.max), true.acf[2,1,]/true.acf[2,1,lag.max+1], col = "red")
plot(ccf.list[[2]], main = "R-CCF between 1 and 2", ylim = range(ccf.list[[2]]$acf))

lines(seq(-lag.max, lag.max), true.acf[2,1,]/true.acf[2,1,lag.max+1], col = "red")
dev.off()

###################################
nsim <- 5e4
###################################

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m


pdf(file = paste("Out/scatter_plot_", nsim, ".pdf", sep = ""), title = paste("Scatter plot, nsim = ", nsim))
plot(mc.chain.list[[1]], xlim = c(-(m/2)*sqrt(diag(target)[1]), (m/2)*sqrt(diag(target)[1])), ylim = c(-(m/2)*sqrt(diag(target)[2]), (m/2)*sqrt(diag(target)[2])), xlab = "X component", ylab = "Y component")
points(mc.chain.list[[2]], col = 2)
points(mc.chain.list[[3]], col = 3)
points(mc.chain.list[[4]], col = 4)
points(mc.chain.list[[5]], col = 5)
dev.off()


acf.list <- combined_acf(mc.chain.list, chain = 1, component = seq(1,p), lag.max = 150, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, chain = 1, component = c(1,2), lag.max = 150, type = "correlation")

pdf(file = paste("Out/acf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim), height = 4)
par(mfrow = c(2,2))

plot(acf.list[[1]][[1]], main = paste("ACF for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[1]][[2]], main = paste("R-ACF for component - 1"), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(acf.list[[2]][[1]], main = paste("ACF for component - 2"), xlim = range(0, lag.max),ylim =  c(0,1))#c(min(min(acf.list[[2]][[1]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[1]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")
plot(acf.list[[2]][[2]], main = paste("R-ACF for component - 2"), xlim = range(0, lag.max), ylim =  c(0,1)) #c(min(min(acf.list[[2]][[2]]$acf), min(true.acf[2,2,])), max(max(acf.list[[2]][[2]]$acf), max(true.acf[2,2,]))))
lines(seq(-lag.max, lag.max), true.acf[2,2,]/true.acf[2,2,lag.max + 1], col = "red")
dev.off()

pdf(file = paste("Out/ccf, n = ", nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered CCF plot, nsim = ", nsim), height = 3)
par(mfrow = c(1,2))
plot(ccf.list[[1]], main = "CCF between 1 and 2", xlim = c(-lag.max, lag.max), ylim = range(ccf.list[[1]]$acf))
lines(seq(-lag.max, lag.max), true.acf[2,1,]/true.acf[2,1,lag.max+1], col = "red")
plot(ccf.list[[2]], main = "R-CCF between 1 and 2", ylim = range(ccf.list[[2]]$acf))
lines(seq(-lag.max, lag.max), true.acf[2,1,]/true.acf[2,1,lag.max+1], col = "red")
dev.off()
