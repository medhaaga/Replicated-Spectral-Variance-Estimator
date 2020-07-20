set.seed(1)
library(fields)
library(graphics)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
library("RColorBrewer")
sourceCpp("lag.cpp")
source("functions.R")


###############################################
##### Model parameters####################
#########################################

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

############3Simuation settings###########

truth <- true.sigma(phi, var = target)
check.pts <- c(1e3, 5e3, 1e4, 5e4, 1e5)
rep <- 50
min <- 5e2
max <- 1e5
step <- 500
conv.pts <- seq(min, max, step)

start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in floor(m/2):1){
  start[i,] <- 2*i*sqrt(diag(target))
  start[m-i+1,] <- -2*i*sqrt(diag(target))
}

###########################################################
############## log Frobenius norm running plots ###########
###########################################################


load(file = paste(paste("Out/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/run_plots", m, sep = "_"), ".pdf", sep = ""), height = 3)
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "coral", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance", 
     ylim = range(c(log(apply(rsv.samp, 3, norm, type = "F"))), log(norm(truth, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="royalblue", lwd = 2)
abline(h = log(norm(truth, type = "F")), col = "green")
legend("bottomright", legend=c("ASV", "RSV"),col=c("coral", "royalblue"), lty=1, lwd = 2)
dev.off()

###########################################################
######### Effective Sample Size running plot##############
##########################################################

load(file = paste(paste("Out/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

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

pdf(file = ("Out/run_plot_ess.pdf"))
plot(conv.pts, mean.asv, type = "l", main = "", xlab = "Simulation size", col = "coral",
     ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)), lwd = 2)
segments(x0 = conv.pts, y0 = lower.asv, x1 = conv.pts, y1 = upper.asv, col = "coral")
lines(conv.pts, mean.rsv, type = "l", col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = lower.rsv, x1 = conv.pts, y1 = upper.rsv, col = "royalblue")
abline(h = (det(target)/det(truth))^(1/p), col = "green3", lwd = 2)
legend("topright", legend=c("ASV", "RSV", "Truth"),col=c("coral", "royalblue", "green3"), lty=1, cex=1, lwd = 2)
dev.off()

##########################################
########### Scatter plots ################
##########################################

master.chain <- array(0, dim = c(5e4,p,m))

for (k in 1:m){
  master.chain[,,k] <- markov.chain(phi, omega, 5e4, start[k,])
}

pdf(file = "Out/scatter_plot.pdf", height = 4)

par(mfrow = c(1,2))

nsim <- 1e3
plot(master.chain[1:nsim,,1], col = 1, xlim = range(master.chain[,1,]), 
     ylim = range(master.chain[,2,]), xlab = "Component-1", ylab = "Component-2", main = "")
for (k in 2:m) {points(master.chain[1:nsim,,k], col = k)}

nsim <- 5e4
plot(master.chain[1:nsim,,1], col = 1, xlim = range(master.chain[,1,]), 
     ylim = range(master.chain[,2,]), xlab = "Component-1", ylab = "Component-2", main = "")
for (k in 2:m) {points(master.chain[1:nsim,,k], col = k)}
dev.off()


######################################################################
############log-Frobenius norm density plots and coverage prob########
#########################################################################


for (j in 1:length(check.pts)){

  nsim <- check.pts[j]
  load(file = paste(paste("Out/out", m, nsim, sep = "_"), ".Rdata", sep = ""))
  
  frob.asv <- density(log(apply(asv.samp, 3, norm, type = "F")))
  frob.rsv <- density(log(apply(rsv.samp, 3, norm, type = "F")))
  
  pdf(file = paste(paste("Out/densities_n", nsim, sep = "_"), ".pdf", sep = ""))
  plot(frob.asv, col = "coral", main = "", xlab = "log-Frobenius norm", ylab = "Density", lwd = 2,
       ylim = c(0, max(frob.asv$y, frob.rsv$y)), xlim = range(frob.asv$x, frob.rsv$x, log(norm(truth, "F"))))
  lines(frob.rsv, col="royalblue", lwd = 2)
  abline(v = log(norm(truth, type = "F")), col = "green3", lwd = 2)
  legend("topleft", legend=c("ASV", "RSV", "Truth"), col=c("coral", "royalblue", "green3"), lty=1, lwd = 2)

  print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
  
  dev.off()
}
