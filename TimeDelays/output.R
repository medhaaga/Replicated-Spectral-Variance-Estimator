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
library(timedelay)
sourceCpp("lag.cpp")
source("functions.R")
load("dataset")

###################################################
### Data and model parameters
lcA <- data[, 1 : 3]
lcB <- data[, c(1, 4, 5)]
m <- 4
delta.start <- c(300, 500, 700, 900)
micro <- 0

###############################################


check.pts <- c(5e2, 1e3, 5e3, 1e4)
r <- length(check.pts)
freq <- 1e2
c.prob <- .95
min <- 5e2
max <- 1e5
step <- 500
conv.pts <- seq(min, max, step)

## 1.) Running plots

load(file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/run_plots_time_delay", sep = "_"), ".pdf", sep = ""))
plot(conv.pts,asv.samp[1,1,], type = "l", col = "red", main = "Variance of time delay", xlab = "Simulation size", ylab = "Variance", ylim = range(rsv.samp[1,1,], asv.samp[1,1,]))
lines(conv.pts, rsv.samp[1,1,], col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
dev.off()

#### 1c.) Determinant
pdf(file = paste(paste("Out/run_plots_det", sep = "_"), ".pdf", sep = ""), height = 3)

det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))

for (i in 1:length(conv.pts)){
  det.rsv[i] <- (det(rsv.samp[,,i]))
  det.asv[i] <- (det(asv.samp[,,i]))
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
#lower.asv <- rep(0, length(conv.pts))
#lower.rsv <- rep(0, length(conv.pts))
#upper.asv <- rep(0, length(conv.pts))
#upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .9)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
  mean.asv[i] <- log(asv[1])
  mean.rsv[i] <- log(rsv[1])
 # lower.asv[i] <- log(asv[2])
#  lower.rsv[i] <- log(rsv[2])
 # upper.asv[i] <- log(asv[3])
  #upper.rsv[i] <- log(rsv[3])
}
plot(conv.pts, mean.asv, type = "l", col = "red", main = "ESS running plot", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(mean.asv, mean.rsv))
lines(conv.pts, mean.rsv, type = "l", col = "blue")
#segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "red")
#segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()



## 3.) Density plots for determinant of ASV and RSV.

pdf(file = "Out/densities.pdf", height = 5)

par(mfrow = c(2,3))
for (j in 1:5){
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
  
  #print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  #print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}
dev.off()

