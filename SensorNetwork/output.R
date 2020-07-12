set.seed(1)
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

###################################################
### Data and model parameters
######################################################################

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors.
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)

# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546,
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0,
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)


m <- 2
start <- rbind(runif(n=8, min=-0.3, max=0), runif(n=8, min=0.7, max=1))
aux <- rbind(runif(n=8, min=-0.3, max=0), runif(n=8, min=0.7, max=1))
j.scale <- rep(1,.08, 4)
truth <- c(0.5748, 0.9069, 0.0991, 0.3651, 0.2578, 0.1350, 0.8546, 0.0392)


###############################################


check.pts <- c(5e3, 1e4, 5e4, 1e5)
r <- length(check.pts)
freq <- 1e2
c.prob <- .95
min <- 1e3
max <- 1e5
step <- 100
conv.pts <- seq(min, max, step)

## 1.) Running plots

load(file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

for (k in 1:p){
  pdf(file = paste("Out/run_plots-", k, ".pdf", sep = ""))
  plot(conv.pts,asv.samp[k,k,], type = "l", col = "red", main = paste("Variance of component -", k), xlab = "Simulation size", ylab = "Variance", ylim = range(rsv.samp[k,k,], asv.samp[k,k,]))
  lines(conv.pts, rsv.samp[k,k,], col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
  dev.off()
}
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

par(mfrow = c(k,k))
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
  
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}
dev.off()

