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

#parameters and initiate values

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
truth <- true.sigma(phi, var = sigma)
true.det <- det(truth)
check.pts <- c(3e2, 5e2, 1e3, 5e3, 1e4)
freq <- 1e3  #100 for now, will change later
rep <- 10
c.prob <- .95
min <- 5e2
max <- 1e4
conv.pts <- seq(min, max, 100)
r <- length(check.pts)
m <- 2
start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in floor(m/2):1){
  start[i,] <- 2*i*sqrt(diag(sigma))
  start[m-i+1,] <- -2*i*sqrt(diag(sigma))
}

#################################
### Bartlett
################################


## 1.) Running plots


load(file = paste(paste("Out/bartlett/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

#Calculate determinants

det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))

for (i in 1:length(conv.pts)){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}

#### 1a.) 5 - components

pdf(file = paste(paste("Out/bartlett/run_plots", m, sep = "_"), ".pdf", sep = ""), height = 3)
par(mfrow = c(2,3))

for (k in 1:p){
  plot(conv.pts,asv.samp[k,k,], type = "l", col = "red", main = paste("Variance of component -  ", k), xlab = "Simulation size", ylab = paste("Sigma_", k, sep = "") , ylim = c(min(rsv.samp[k,k,]), max(truth[k,k], max(rsv.samp[k,k,]))))
  lines(conv.pts, rsv.samp[k,k,], col="blue")
  abline(h = truth[k,k], col = "green")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
}

#### 1c.) Determinant

plot(conv.pts,det.asv, type = "l", col="red", main = "Running plot of determinant", xlab = "Simulation size", ylab = "Determinant" ,ylim = c(min(det.rsv), max(det(truth), max(det.rsv))))
lines(conv.pts,det.rsv, col="blue")
abline(h = det(truth), col = "green")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
dev.off()

#### 1d.) ESS


mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
#lower.asv <- rep(0, length(conv.pts))
#lower.rsv <- rep(0, length(conv.pts))
#upper.asv <- rep(0, length(conv.pts))
#upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .9)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  #lower.asv[i] <- asv[2]
  #lower.rsv[i] <- rsv[2]
  #upper.asv[i] <- asv[3]
  #upper.rsv[i] <- rsv[3]
}

pdf(file = paste(paste("Out/bartlett/run_plot_ess", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts, mean.asv, type = "l", col = "red", main = "Running plot for ESS/mn", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts, mean.rsv, type = "l", col = "blue")
abline(h = (det(sigma)/det(truth))^(1/p), col = "green")
legend("topright", legend=c("ASV", "RSV", "Truth"),col=c("red", "blue", "green"), lty=1, cex=1)
dev.off()


## 2.) Scatter plots and trace plots

master.chain <- array(0, dim = c(max(check.pts),p,m))

for (k in 1:m){
  master.chain[,,k] <- markov.chain(phi, omega, max(check.pts), start[k,])
}

#### 2a.) Scatter plots
pdf(file = paste(paste("Out/bartlett/scatter_plot", m, sep = "_"), ".pdf", sep = ""))
par(mfrow = c(2,2))
for (j in 1:4){

  nsim <- check.pts[j]
  chain <- master.chain[1:nsim,,]
  plot(chain[,1,1], chain[,5,1], col = 2, xlim = c(-m*sqrt(diag(sigma))[1], m*sqrt(diag(sigma))[1]), ylim = c(-m*sqrt(diag(sigma))[1], m*sqrt(diag(sigma))[1]), main = paste("Sample size = ", nsim), xlab = "Component - 1", ylab = "Component - 5")
  for (k in 2:m){
    points(chain[,1,k], chain[,p,k], col = k+1)
  }
}
dev.off()

#### 2b.) Trace plots

for (j in 1:r){
  nsim = check.pts[j]
  chain <- master.chain[1:nsim,,]
  pdf(file = paste(paste("Out/bartlett/trace", m, nsim, sep = "_"), ".pdf", sep = ""))
  par(mfcol = c(p,m))
  for(k in 1:m){
    for (t in 1:p){
      plot.ts(chain[,t,k], main = paste("Component - ", t))
    }
  }
  dev.off()
}

## 3.) Density plots for determinant of ASV and RSV

pdf(file = paste(paste("Out/bartlett/densities",m, sep = "_"), ".pdf", sep = ""))
par(mfrow = c(ceil(r/2),2))
for (j in 1:r){  

  #Comparitative density plots of ASV and RSV
  nsim <- check.pts[j]
  load(file = paste(paste("Out/bartlett/out", m, nsim, sep = "_"), ".Rdata", sep = ""))

  #Calculate determinants
  det.rsv <- rep(0,freq)
  det.asv <- rep(0,freq)
  for (f in 1:freq){
    det.rsv[f] <- det(rsv.samp[,,f])
    det.asv[f] <- det(asv.samp[,,f])
  }

  plot(density(det.asv), col="red", main = paste("Determinant, nsim = ", nsim), xlab = "Determinant", ylab = "Density", xlim = c(min(det.rsv), max(det(truth), max(det.rsv))))
  lines(density(det.rsv), col="blue")
  abline(v = det(truth), col = "green")
  legend("topright", legend=c("ASV", "RSV", "Truth"),col=c("red", "blue", "green"), lty=1, cex=1.2)


  print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

dev.off()

#################################
### tukey
################################


## 1.) Running plots


load(file = paste(paste("Out/tukey/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

#Calculate determinants

det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))

for (i in 1:length(conv.pts)){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}


pdf(file = paste(paste("Out/tukey/run_plots", m, sep = "_"), ".pdf", sep = ""), height = 3)
par(mfrow = c(2,3))

#### 1a.) 5 - components

for (k in 1:p){
  plot(conv.pts,asv.samp[k,k,], type = "l", col = "red", main = paste("Variance of component -  ", k), xlab = "Simulation size", ylab = paste("Sigma_", k, sep = "") , ylim = c(min(rsv.samp[k,k,]), max(truth[k,k], max(rsv.samp[k,k,]))))
  lines(conv.pts, rsv.samp[k,k,], col="blue")
  abline(h = truth[k,k], col = "green")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
}

#### 1c.) Determinant

plot(conv.pts,det.asv, type = "l", col="red", main = "Running plot of determinant", xlab = "Simulation size", ylab = "Determinant" ,ylim = c(min(det.rsv), max(det(truth), max(det.rsv))))
lines(conv.pts,det.rsv, col="blue")
abline(h = det(truth), col = "green")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
dev.off()

#### 1d.) ESS


mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
#lower.asv <- rep(0, length(conv.pts))
#lower.rsv <- rep(0, length(conv.pts))
#upper.asv <- rep(0, length(conv.pts))
#upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .9)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  #lower.asv[i] <- asv[2]
  #lower.rsv[i] <- rsv[2]
  #upper.asv[i] <- asv[3]
  #upper.rsv[i] <- rsv[3]
}

pdf(file = paste(paste("Out/tukey/run_plot_ess", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts, mean.asv, type = "l", col = "red", main = "Running plot for ESS/mn", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts, mean.rsv, type = "l", col = "blue")
abline(h = (det(sigma)/det(truth))^(1/p), col = "green")
legend("topright", legend=c("ASV", "RSV", "Truth"),col=c("red", "blue", "green"), lty=1, cex=1)
dev.off()


## 2.) Scatter plots and trace plots

master.chain <- array(0, dim = c(max(check.pts),p,m))

for (k in 1:m){
  master.chain[,,k] <- markov.chain(phi, omega, max(check.pts), start[k,])
}

#### 2a.) Scatter plots
pdf(file = paste(paste("Out/tukey/scatter_plot", m, sep = "_"), ".pdf", sep = ""))
par(mfrow = c(2,2))
for (j in 1:4){
  
  nsim <- check.pts[j]
  chain <- master.chain[1:nsim,,]
  plot(chain[,1,1], chain[,5,1], col = 2, xlim = c(-m*sqrt(diag(sigma))[1], m*sqrt(diag(sigma))[1]), ylim = c(-m*sqrt(diag(sigma))[1], m*sqrt(diag(sigma))[1]), main = paste("Sample size = ", nsim), xlab = "Component - 1", ylab = "Component - 5")
  for (k in 2:m){
    points(chain[,1,k], chain[,p,k], col = k+1)
  }
}
dev.off()

#### 2b.) Trace plots

for (j in 1:r){
  nsim = check.pts[j]
  chain <- master.chain[1:nsim,,]
  pdf(file = paste(paste("Out/tukey/trace", m, nsim, sep = "_"), ".pdf", sep = ""))
  par(mfcol = c(p,m))
  for(k in 1:m){
    for (t in 1:p){
      plot.ts(chain[,t,k], main = paste("Component - ", t))
    }
  }
  dev.off()
}

## 3.) Density plots for determinant of ASV and RSV

pdf(file = paste(paste("Out/tukey/densities",m, sep = "_"), ".pdf", sep = ""))
par(mfrow = c(ceil(r/2),2))
for (j in 1:r){  
  
  #Comparitative density plots of ASV and RSV
  nsim <- check.pts[j]
  load(file = paste(paste("Out/tukey/out", m, nsim, sep = "_"), ".Rdata", sep = ""))
  
  #Calculate determinants
  det.rsv <- rep(0,freq)
  det.asv <- rep(0,freq)
  for (f in 1:freq){
    det.rsv[f] <- det(rsv.samp[,,f])
    det.asv[f] <- det(asv.samp[,,f])
  }
  
  plot(density(det.asv), col="red", main = paste("Determinant, nsim = ", nsim), xlab = "Determinant", ylab = "Density", xlim = c(min(det.rsv), max(det(truth), max(det.rsv))))
  lines(density(det.rsv), col="blue")
  abline(v = det(truth), col = "green")
  legend("topright", legend=c("ASV", "RSV", "Truth"),col=c("red", "blue", "green"), lty=1, cex=1.2)
  
  
  print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

dev.off()