set.seed(1)
library(fields)
library(scatterplot3d)
library(cubature)
library(graphics)
library(pracma)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
sourceCpp("lag.cpp")
source("functions.R")
library("MCMCpack")
load("mida.rda")
library(rep.acf.ccf)


###############################################

c0 <- 13
d0 <- 1
bars <- 6
p <- bars+1
m <- 2
check.pts <- c(5e2, 1e3, 2e3, 5e3, 1e4, 2e4)
r <- length(check.pts)
freq <- 1e2
c.prob <- .95
min <- 5e2
max <- 1e5
step <- 500
conv.pts <- seq(min, max, step)

## 1.) Running plots

load(file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/run_plots_components", sep = "_"), ".pdf", sep = ""), height = 3)
par(mfrow=c(4,2))

for (i in 1:(bars+1)){

  plot(conv.pts,asv.samp[i,i,], type = "l", col = "red", main = paste("Variance of component -", i), xlab = "Simulation size", ylab = "Variance", ylim = range(rsv.samp[i,i,]))
  lines(conv.pts, rsv.samp[i,i,], col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)

}


#### 1c.) Determinant
pdf(file = paste(paste("Out/run_plots_det", sep = "_"), ".pdf", sep = ""), height = 3)

det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))

for (i in 1:length(conv.pts)){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}

plot(conv.pts,det.asv, type = "l", col="red", main = paste("Determinant"), xlab = "Simulation size", ylab = "Determinant", ylim = range(det.rsv))
lines(conv.pts,det.rsv, col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
dev.off()

#### 1d.) ESS

pdf(file = paste("Out/run_plot_ess.pdf", sep = "_"), height = 4)

par(mfrow = c(1,1))
mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .9)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  lower.asv[i] <- asv[2]
  lower.rsv[i] <- rsv[2]
  upper.asv[i] <- asv[3]
  upper.rsv[i] <- rsv[3]
}
plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "red", main = "ESS running plot", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "blue")
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "red")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)



## 2.) ACF and trace plots

chain <- array(0, dim = c(1e5, p, m))
for (j in 1:m){
  chain[,,j] <- matrix(MCMCpoissonChange(mida~ 1, m = bars, c0=c0, d0=d0, marginal.likelihood="none", mcmc = 1e5, burnin = 0, thin = 1, seed = (499*j), nrow = 1e5, ncol = p))
}

#### 2a.) Trace plots
for(i in 1:p){
  pdf(file = paste(paste("Out/trace_component", i, sep = "_"), ".pdf", sep = ""))
  par(mfrow = c(4,2))
  for (j in 3:r){
    nsim = check.pts[j]
    plot.ts(chain[1:nsim,i,1], main = paste("Chain - 1, nsim -", nsim))
    plot.ts(chain[1:nsim,i,2], main = paste("Chain - 2, nsim -", nsim))
  }
  dev.off()
}

#### 2b.) ACF-RACF comparison for component-4 of chain-1

ncrop <- 1e3
foo1 <- list(chain[1:ncrop,,1], chain[1:ncrop,,2])
racf_crop <- combined_acf(foo1, chain = 1, component = c(3), lag.max = 100, type = "correlation")
foo2 <- list(chain[,,1], chain[,,2])
racf <- combined_acf(foo2, chain = 1, component = c(3), lag.max = 100, type = "correlation")

pdf(file = "Out/acf_racf.pdf")
par(mfrow = c(2,2))
plot(racf_crop[[1]][[1]], main = "ACF for nsim = 1000")
plot(racf_crop[[1]][[2]], main = "R-ACF for nsim = 1000")
plot(racf[[1]][[1]], main = "ACF for nsim = 100000")
plot(racf[[1]][[2]], main = "R-ACF for nsim = 100000")
dev.off()


## 3.) Density plots for determinant of ASV and RSV.

pdf(file = "Out/densities.pdf", height = 5)

par(mfrow = c(2,3))
for (j in 1:3){
  nsim <- check.pts[j]
  load(file = paste(paste("Out/out", nsim, sep = "_"), ".Rdata", sep = ""))

  det.rsv <- rep(0,freq)
  det.asv <- rep(0,freq)
  for (k in 1:freq){
    det.rsv[k] <- det(rsv.samp[,,k])
    det.asv[k] <- det(asv.samp[,,k])
  }
  plot(density(det.asv), col="red", main = paste(" n =", nsim), xlab = "Determinant")
  lines(density(det.rsv), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}
dev.off()

