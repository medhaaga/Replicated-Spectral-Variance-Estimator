set.seed(1)
library(multichainACF)


log.density <- function(x, p, mu1, mu2, sd1, sd2){
  return(log(p*dnorm(x, mean = mu1, sd = sd1) + (1-p)*dnorm(x, mean = mu2, sd= sd2)))
}


mh.mcmc <- function(start, p, mu1, mu2, sd1, sd2, N, h){
  X <- rep(0, N)
  X[1] <- start
  for (i in 2:N){
    prop <- rnorm(1, mean = X[i-1], sd = 2)
    ratio <- log.density(prop, p, mu1, mu2, sd1, sd2) - log.density(X[i-1], p, mu1, mu2, sd1, sd2)
    if(runif(1) < exp(ratio))
      X[i] <- prop else
        X[i] <- X[i-1]
  }
  return (X)
}


N <- 1e5
p <- 0.7
mu1 <- -5
mu2 <- 5
sd1 <- 1
sd2 <- 0.5


chain1 <- mh.mcmc(start = -3, p, mu1, mu2, sd1, sd2, N, 1)
chain2 <- mh.mcmc(start = 3, p, mu1, mu2, sd1, sd2, N, 1)
mean(chain1)
mean(chain2)
mc.chain.list <- list(as.matrix(chain1), as.matrix(chain2))
save(mc.chain.list, file = "gaussian-two_chains.Rdata")

load(file = "gaussian-two_chains.Rdata")


##################################
######Trace plots for n=1e4#######
##################################

x <- seq(-10, 10, length = 1e3)
pdf(file = "TargetTrace_n1e4.pdf", height = 5, width = 6)
plot(x, 20000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x",cex.lab=1.2, ylab = "", ylim = c(-10000,6000), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1[1:1e4], y = seq(-1, -1e4, -1), col = "lightskyblue", xlab = "", ylab = "", type = "l", ylim = c(-1e4, 6e3), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2[1:1e4], y = seq(-1, -1e4, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1, cex=1.2)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=.7, lwd=2)
dev.off()

##################################
######Trace plots for n=1e5 ######
##################################


pdf(file = "TargetTrace_n1e5.pdf", height = 5, width = 6)
plot(x, 200000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x", cex.lab=1.2, ylab = "", ylim = c(-1e5,6e4), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1, y = seq(-1, -1e5, -1), col = "lightskyblue", xlab = "", ylab = "", cex.lab=1.2, type = "l", ylim = c(-1e5, 6e4), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2, y = seq(-1, -1e5, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1, cex=1.2)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.7, lwd=2)
dev.off()


##############################################
########## ACF and G-ACF ####################
#############################################

m <- 2
lag.max <- 50
nsim1 <- 1e4
nsim2 <- 1e5
########################################

x <- list()
y <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:nsim1,])
  y[[i]] <- as.matrix(mc.chain.list[[i]][1:nsim2,])
}

local.acf1 <- globalACF(x, chains = c(1), component = 1, mean = "local", lag.max = lag.max, type = "correlation", avg = FALSE, plot = FALSE)[[1]]
local.acf2 <- globalACF(y, chains = c(1), component = 1, mean = "local", lag.max = lag.max, type = "correlation", avg = FALSE, plot = FALSE)[[1]]
global.acf1 <- globalACF(x, chains = c(1), component = 1, lag.max = lag.max, type = "correlation", avg = FALSE, plot = FALSE)[[1]]
global.acf2 <- globalACF(y, chains = c(1), component = 1, lag.max = lag.max, type = "correlation", avg = FALSE, plot = FALSE)[[1]]

pdf(file = "gaussian-acf_hist.pdf", width = 10, height= 4)
par(mfrow = c(1,2))
plot(as.matrix(local.acf1$acf), type = 'h', ylab = "Autocorrelation", xlab = "Lag", cex.lab=1.2)
lines(as.matrix(local.acf2$acf), type = 'l', col = "steelblue1", lwd = 2)
plot(as.matrix(global.acf1$acf), type = 'h', ylim = c(min(local.acf1$acf), 1), ylab = "Autocorrelation", xlab = "Lag", cex.lab=1.2)
lines(as.matrix(global.acf2$acf), type = 'l', col = "steelblue1", lwd = 2)
dev.off()
