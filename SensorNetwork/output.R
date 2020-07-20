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


m <- 5
start <- rbind(runif(n=8, min=-0.3, max=0), runif(n=8, min=0.7, max=1))
aux <- rbind(runif(n=8, min=-0.3, max=0), runif(n=8, min=0.7, max=1))
j.scale <- rep(1,.08, 4)
truth <- c(0.5748, 0.9069, 0.0991, 0.3651, 0.2578, 0.1350, 0.8546, 0.0392)


###############################################


check.pts <- c(5e3, 1e4, 5e4, 1e5)
r <- length(check.pts)
freq <- 1e2
c.prob <- .95
min <- 500
max <- 2e5
step <- 500
conv.pts <- seq(min, max, step)

## 1.) Running plots
p <- 8
load(file = paste("Out/conv_data_m", m, "_min", min, "_max", max, ".Rdata", sep = ""))

# for (k in 1:p){
#   a <- lapply(asv, function(x) x[k,k, ])
#   r <- lapply(rsv, function(x) x[k,k, ])
#   a <- Reduce("+", a)/length(a)
#   r <- Reduce("+", r)/length(r)
#   pdf(file = paste("Out/run_plots-", k, ".pdf", sep = ""))
#   plot(conv.pts, a, type = "l", col = "red", main = paste("Variance of component -", k), xlab = "Simulation size", ylab = "Variance", ylim = range(a, r))
#   lines(conv.pts, r, col="blue")
#   legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
#   dev.off()
# }


ind <- seq(1,length(conv.pts), by = 2)

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)
  pdf(file = paste("Out/run_plots-Frob.pdf", sep = ""), height = 5, width = 5)
  plot(conv.pts, a, type = "l", col = "red", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)) )
  lines(conv.pts, r, col="blue")
  segments(x0 = conv.pts[ind], y0 = (a - se.a)[ind], y1 = (a + se.a)[ind], col = adjustcolor("red", alpha.f = .70))
segments(x0 = conv.pts[ind], y0 = (r - se.r)[ind], y1 = (r + se.r)[ind], col = adjustcolor("blue", alpha.f = .70))

  legend("bottomright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1)
  dev.off()

#### 1c.) Determinant
# pdf(file = paste(paste("Out/run_plots_det", sep = "_"), ".pdf", sep = ""), height = 3)

# det.rsv <- rep(0,length(conv.pts))
# det.asv <- rep(0,length(conv.pts))

#   a <- lapply(asv, function(x) det(x[, ,i])^(1/p))
#   r <- lapply(rsv, function(x) det(x[, ,i])^(1/p))
# for (i in 1:length(conv.pts)){

#   a <- lapply(asv, function(x) det(x[, ,i])^(1/p))
#   r <- lapply(rsv, function(x) det(x[, ,i])^(1/p))
#   # a <- Reduce("+", a)/length(a)
#   # r <- Reduce("+", r)/length(r)
#   det.asv[i] <- Reduce("+", a)/length(a)
#   det.rsv[i] <- Reduce("+", r)/length(r)
# }

# plot(conv.pts,det.asv, type = "l", col="red", main = paste("Determinant"), xlab = "Simulation size", ylab = "Determinant", ylim = range(c(det.asv, det.rsv)))
# lines(conv.pts,det.rsv, col="blue")

# legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
# dev.off()

#### 1d.) ESS

pdf(file = paste("Out/run_plot_ess.pdf", sep = "_"), height = 5, width = 5)

# par(mfrow = c(k,k))
# mean.asv <- rep(0, length(conv.pts))
# mean.rsv <- rep(0, length(conv.pts))
#lower.asv <- rep(0, length(conv.pts))
#lower.rsv <- rep(0, length(conv.pts))
#upper.asv <- rep(0, length(conv.pts))
#upper.rsv <- rep(0, length(conv.pts))
# for (i in 1:length(conv.pts)){
#   asv <- confidence_interval(ess.asv.samp[[i]], .9)
#   rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
#   mean.asv[i] <- log(asv[1])
#   mean.rsv[i] <- log(rsv[1])
  # lower.asv[i] <- log(asv[2])
  #  lower.rsv[i] <- log(rsv[2])
  # upper.asv[i] <- log(asv[3])
  #upper.rsv[i] <- log(rsv[3])
# }
# plot(conv.pts, mean.asv, type = "l", col = "red", main = "ESS running plot", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(mean.asv, mean.rsv))
# lines(conv.pts, mean.rsv, type = "l", col = "blue")
#segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "red")
#segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")

  a <- lapply(ess.asv, log)
  r <- lapply(ess.rsv, log)
  se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
  se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
  a <- Reduce("+", a)/length(ess.asv)
  r <- Reduce("+", r)/length(ess.rsv)

plot(conv.pts, a, type = "l", col = "red", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r))
lines(conv.pts, r, col = "blue")
  segments(x0 = conv.pts[ind], y0 = (a - se.a)[ind], y1 = (a + se.a)[ind], col = adjustcolor("red", alpha.f = .70))
segments(x0 = conv.pts[ind], y0 = (r - se.r)[ind], y1 = (r + se.r)[ind], col = adjustcolor("blue", alpha.f = .70))

legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()



## 3.) Density plots for determinant of ASV and RSV.

# pdf(file = "Out/densities.pdf", height = 5)

# par(mfrow = c(2,3))
# for (j in 1:5){
#   nsim <- check.pts[j]
#   load(file = paste(paste("Out/out", nsim, sep = "_"), ".Rdata", sep = ""))
  
#   det.rsv <- rep(0,freq)
#   det.asv <- rep(0,freq)
#   for (k in 1:freq){
#     det.rsv[k] <- det(rsv.samp[,,k])
#     det.asv[k] <- det(asv.samp[,,k])
#   }
#   plot(density(det.asv), col="red", main = paste(" n =", nsim), xlab = "Determinant")
#   lines(density(det.rsv), col="blue")
#   legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  
#   print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
#   print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
# }
# dev.off()

