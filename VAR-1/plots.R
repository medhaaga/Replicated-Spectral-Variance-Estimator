set.seed(1)
library(fields)
library(scatterplot3d)
library(graphics)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
sourceCpp("lag.cpp")
source("functions.R")


###############################################
phi <- diag(c(.99, .01))
omega <- diag(2)
m = 2
p <- ncol(phi)

#sims for plotting densities and calculating coverage

check.pts <- c(3e2, 5e2, 1e3, 5e3, 1e4)
freq <- 1e2  #100 for now, will change later
rep <- 10
c.prob <- .95
min <- 5e2
max <- 1e4
conv.pts <- seq(min, max, 100)
r <- length(check.pts)
start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in floor(m/2):1){
  start[i,] <- 2*i*sqrt(diag(truth))
  start[m-i+1,] <- -2*i*sqrt(diag(truth))
}


## 1.) Running plots
load(file = paste(paste("Out/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

#Calculate determinants

det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))

for (i in 1:length(conv.pts)){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}

#### 1a.) Sigma_11
pdf(file = paste(paste("Out/run_plots", m, sep = "_"), ".pdf", sep = ""), height = 3)

par(mfrow = c(1,3))
plot(conv.pts,asv.samp[1,1,], type = "l", col = "red", main = "Running plot of Sigma_11", xlab = "Simulation size", ylab = "Sigma_11", ylim = range(rsv.samp[1,1,]))
lines(conv.pts, rsv.samp[1,1,], col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)


#### 1b.) Sigma_22

plot(conv.pts,asv.samp[2,2,], type = "l", col = "red", main = "Running plot of Sigma_22", xlab = "Simulation size", ylab = "Sigma_22", ylim = range(rsv.samp[2,2,]))
lines(conv.pts, rsv.samp[2,2,], col="blue")


#### 1c.) Determinant

plot(conv.pts,det.asv, type = "l", col="red", main = "Running plot of determinant", xlab = "Simulation size", ylab = "Determinant" ,ylim = range(det.rsv))
lines(conv.pts,det.rsv, col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()

#### 1d.) ESS

ess.asv <- rep(0, length(conv.pts))
ess.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  ess.rsv[i] <- ess.rsv.samp[[i]][1]
  ess.asv[i] <- ess.asv.samp[[i]][1]
  lower.asv[i] <- ess.asv.samp[[i]][2]
  lower.rsv[i] <- ess.rsv.samp[[i]][2]
  upper.asv[i] <- ess.asv.samp[[i]][3]
  upper.rsv[i] <- ess.rsv.samp[[i]][3]
}


pdf(file = paste(paste("Out/run_plot_ess", m, sep = "_"), ".pdf", sep = ""))
plotCI(conv.pts, ess.asv, li = lower.asv, ui = upper.asv, col = "red")
par(new=T)
plotCI(conv.pts, ess.rsv, li = lower.rsv, ui = upper.rsv, col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()


## 2.) Density plots for Sigma11, Sigma22 and determinant.
##     Loop over check points for different values of nsim

for (j in 1:r){

  #### 2a.) Scatter plot of m Markov chains
  nsim <- check.pts[j]
  chain <- array(0, dim = c(nsim,p,m))

  for (k in 1:m){
    chain[,,k] <- markov.chain(phi, omega, nsim, start[k,])
  }

  pdf(file = paste(paste("Out/scatter_plot", m, nsim, sep = "_"), ".pdf", sep = ""))
  plot(chain[,1,1], chain[,2,1], col = 2, xlim = c(-m*sqrt(diag(truth))[1], m*sqrt(diag(truth))[1]), ylim = c(-m*sqrt(diag(truth))[1], m*sqrt(diag(truth))[1]), main = paste("Scatter plot of m = ", m,  "Markov chains, n = ", nsim), xlab = "X", ylab = "Y")
  for (p in 2:m){
    points(chain[,1,p], chain[,2,p], col = p+1)
  }

  dev.off()

  pdf(file = paste(paste("Out/trace", m, nsim, sep = "_"), ".pdf", sep = ""))
  par(mfrow = c(m,p))
  for(k in 1:m){
    plot.ts(chain[,1,k], main = paste("x component of chain -", k))
    plot.ts(chain[,2,k], main = paste("y component of chain -", k))
  }
  dev.off()

  #Comparitative density plots of ASV and RSV
  load(file = paste(paste("Out/out", m, nsim, sep = "_"), ".Rdata", sep = ""))

  #Calculate determinants
  det.rsv <- rep(0,freq)
  det.asv <- rep(0,freq)
  for (j in 1:freq){
    det.rsv[j] <- det(rsv.samp[,,j])
    det.asv[j] <- det(asv.samp[,,j])
  }

  #### 2b.) Sigma_11

  pdf(file = paste(paste("Out/densities",m, nsim, sep = "_"), ".pdf", sep = ""), height = 3)
  par(mfrow = c(1,3))

  plot(density(asv.samp[1,1,]), col="red", main = paste("Sigma_11, nsim = ", nsim))
  lines(density(rsv.samp[1,1,]), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

  #### 2c.) Sigma_22

  plot(density(asv.samp[2,2,]), col="red", main = paste("Sigma_22, nsim = ", nsim))
  lines(density(rsv.samp[2,2,]), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

  #### 2d.) Determinant

  plot(density(det.asv), col="red", main = paste("Determinant, nsim = ", nsim))
  lines(density(det.rsv), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()

  print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}


