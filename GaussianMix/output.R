set.seed(1)
library(rep.acf.ccf)


log.density <- function(x, p, mu1, mu2, sd1, sd2){
  return(log(p*dnorm(x, mean = mu1, sd = sd1) + (1-p)*dnorm(x, mean = mu2, sd= sd2)))
}


mh.mcmc <- function(start, p, mu1, mu2, sd1, sd2, N, h){

  X <- rep(0, N)
  acc <- 0
  X[1] <- start
  for (i in 2:N){
    prop <- rnorm(1, mean = X[i-1], sd = 2)
    ratio <- log.density(prop, p, mu1, mu2, sd1, sd2) - log.density(X[i-1], p, mu1, mu2, sd1, sd2)

    if(runif(1) < exp(ratio)){
      X[i] <- prop
      acc <- acc+1
    } else {
      X[i] <- X[i-1]
    }
  }
  print(acc/N)
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

##################################
######Trace plots for n=1e4#######
##################################

x <- seq(-10, 10, length = 1e3)
pdf(file = "TargetTrace_n1e4.pdf", height = 5, width = 6)
plot(x, 20000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x", ylab = "", ylim = c(-10000,6000), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1[1:1e4], y = seq(-1, -1e4, -1), col = "lightskyblue", xlab = "", ylab = "", type = "l", ylim = c(-1e4, 6e3), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2[1:1e4], y = seq(-1, -1e4, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.7, lwd=2)
dev.off()

##################################
######Trace plots for n=1e5 ######
##################################


pdf(file = "TargetTrace_n1e5.pdf", height = 5, width = 6)
plot(x, 200000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x", ylab = "", ylim = c(-1e5,6e4), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1, y = seq(-1, -1e5, -1), col = "lightskyblue", xlab = "", ylab = "", type = "l", ylim = c(-1e5, 6e4), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2, y = seq(-1, -1e5, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.7, lwd=2)
dev.off()


###############################################
##############Density plots ####################
###############################################

pdf(file = "density.pdf", height = 5, width = 5)
x <- seq(-10, 10, .01)
plot(density(chain1[1:1e4]), col = "lightskyblue", type = "l", lwd=2, xlim = range(chain1, chain2), ylim = c(0, 1), xlab = "x", ylab = "density", main = "")
polygon(density(chain1[1:1e4]), col = rgb(135, 206, 250, max = 255, alpha = 100, names = "lsb50"))
lines(density(chain2[1:1e4]), col = "plum3", lwd=2)
polygon(density(chain2[1:1e4]), col = rgb(205, 150, 205, max = 255, alpha = 100, names = "p350"))
lines(x, exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.8, lwd=2)

dev.off()


##############################################
########## ACF and G-ACF ####################
#############################################

m <- 2
mc.chain.list <- list(as.matrix(chain1), as.matrix(chain2))

######### ncrop = 1e4###################

####################################
nsim <- 1e4
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:nsim,])

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, 1], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, 1], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = "acf_n1e4.pdf", height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf, type = "h", xlab = "Lag", ylab = "Autocorrelation", main = "", ylim = c(0,1))
plot(global.acf, type = "h", xlab = "Lag", ylab = "Autocorrelation", main = "", ylim = c(0,1))
dev.off()

################## ncrop = 1e5 #########################

####################################
nsim <- 1e5
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:nsim,])

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, 1], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, 1], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = "acf_n1e4.pdf", height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf, type = "h", xlab = "Lag", ylab = "Autocorrelation", main = "", ylim = c(0,1))
plot(global.acf, type = "h", xlab = "Lag", ylab = "Autocorrelation", main = "", ylim = c(0,1))
dev.off()
