###################################################################
##This code produces and saves all the plots in pdf and prints the
##contents of tables using R objects generated in simulations.R
###################################################################

set.seed(1)
library(cubature)
library(graphics)
library(pracma)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
library("RColorBrewer")
sourceCpp("lag.cpp")
source("functions.R")


########################################################
#####Visualization through perspective plots############
########################################################

#Perspective plots for both simulation settings

A <- 1
B <- 3
C <- 8

pdf(file = paste("Out/", A, B, C, "/3d_density_plot.pdf", sep = "_"))
samples <- perspective(A, B, C, 100)
persp(samples$x, samples$y, samples$z, main="", col = "white", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()

A <- 1
B <- 10
C <- 7

pdf(file = paste("Out/", A, B, C, "/3d_density_plot.pdf", sep = "_"))
samples <- perspective(A, B, C, 200)
persp(samples$x, samples$y, samples$z, main="", col = "white", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()


########################################################
########################################################
####### Scatter plots for both simulation settings######
########################################################
########################################################


m <- 2

########### A=1, B=3, C=8##############

set.seed(1)  # use seed for reproducible results
A <- 1
B <- 3
C <- 8
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

chain <- array(0, dim = c(5e4, p, m))
chain[,,1] <- markov.chain(A, B, C, 5e4, start[1,])
chain[,,2] <- markov.chain(A, B, C, 5e4, start[2,])

pdf(file = paste(paste("Out/", A, B, C, "/scatter_plot", m, sep = "_"), ".pdf", sep = ""), height = 4)
par(mfrow = c(1,2))
plot(chain[1:1e3,,1], col = "coral", xlim = range(chain[,1,]), ylim = range(chain[,2,]))
points(chain[1:1e3,,2], col = "blue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("coral", "blue"), lty=1)

plot(chain[,,1], col = "coral", xlim = range(chain[,1,]), ylim = range(chain[,2,]))
points(chain[,,2], col = "blue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("coral", "blue"), lty=1)
dev.off()


########### A=1, B=10, C=7 ##############

set.seed(1)  # use seed for reproducible results
A <- 1
B <- 10
C <- 7
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

chain <- array(0, dim = c(1e4, p, m))
chain[,,1] <- markov.chain(A, B, C, 1e4, start[1,])
chain[,,2] <- markov.chain(A, B, C, 1e4, start[2,])

pdf(file = paste(paste("Out/", A, B, C, "/scatter_plot", m, sep = "_"), ".pdf", sep = ""), height = 4)
par(mfrow = c(1,2))
plot(chain[1:1e3,,1], col = "coral", xlim = range(chain[,1,]), ylim = range(chain[,2,]))
points(chain[1:1e3,,2], col = "blue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("coral", "blue"), lty=1)

plot(chain[,,1], col = "coral", xlim = range(chain[,1,]), ylim = range(chain[,2,]))
points(chain[,,2], col = "blue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("coral", "blue"), lty=1)
dev.off()


###################################################################
###################################################################
########## Running plots for log Frobenius norm ###################
###################################################################
###################################################################

###Running plots for both the simulation settings

##########################################
######### A=1, B=3, C=8 ##################
##########################################

A <- 1
B <- 3
C <- 8

min <- 5e2
max <- 1e5 
step <- 5e2
conv.pts <- seq(min, max, step)

##### two chains ######

m <- 2
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/run_plots", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts[1:100],log(apply(asv.samp, 3, norm, type = "F"))[1:100], type = "l", col = "coral", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance", 
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts[1:100], log(apply(rsv.samp, 3, norm, type = "F"))[1:100], col="blue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, lwd = 2)
dev.off()

##### five chains ######

m <- 5
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/run_plots", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "coral", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance", 
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="blue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, lwd = 2)
dev.off()

###########################################
######### A=1, B=10, C=7 ##################
###########################################

A <- 1
B <- 10
C <- 7

min <- 5e2
max <- 5e4 
step <- 5e2
conv.pts <- seq(min, max, step)

##### two chains ######

m <- 2
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/run_plots", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "coral", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance", 
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="blue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, lwd = 2)
dev.off()

##### five chains ######

m <- 5
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/run_plots", m, sep = "_"), ".pdf", sep = ""))
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "coral", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance", 
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="blue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, lwd = 2)
dev.off()

##########################################################
##########################################################
#######Running plots for Effective Sample Size############
##########################################################
##########################################################

############ESS for both simulation settings##############

###########################################
############# A=1, B=3, C=8 ###############
###########################################

A <- 1
B <- 3 
C <- 8
min <- 5e2
max <- 1e5
step <- 5e2
### convergence is achieved at 50000 iterations only.
conv.pts <- seq(min, 5e4, step)

#### two chains

m <- 2
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste("Out/", A, B, C, "/run_plot_ess_2.pdf", sep = "_"))
mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .95)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .95)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  lower.asv[i] <- asv[2]
  lower.rsv[i] <- rsv[2]
  upper.asv[i] <- asv[3]
  upper.rsv[i] <- rsv[3]
}

plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "coral", main = "Two parallel chains", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "blue")
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "coral")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"), col=c("coral", "blue"), lty=1, cex=1.2)
dev.off()

#### five chains

m <- 5
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste("Out/", A, B, C, "/run_plot_ess_5.pdf", sep = "_"))
mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .95)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .95)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  lower.asv[i] <- asv[2]
  lower.rsv[i] <- rsv[2]
  upper.asv[i] <- asv[3]
  upper.rsv[i] <- rsv[3]
}

plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "coral", main = "Two parallel chains", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "blue")
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "coral")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, cex=1.2)
dev.off()


#############################################
############# A=1, B=10, C=7 ###############
############################################

A <- 1
B <- 10
C <- 7
min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)

#### two chains

m <- 2
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste("Out/", A, B, C, "/run_plot_ess_2.pdf", sep = "_"))
mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .95)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .95)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  lower.asv[i] <- asv[2]
  lower.rsv[i] <- rsv[2]
  upper.asv[i] <- asv[3]
  upper.rsv[i] <- rsv[3]
}

plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "coral", main = "Two parallel chains", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "blue")
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "coral")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, cex=1.2)
dev.off()

#### five chains

m <- 5
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste("Out/", A, B, C, "/run_plot_ess_2.pdf", sep = "_"))
mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))
for (i in 1:length(conv.pts)){
  asv <- confidence_interval(ess.asv.samp[[i]], .95)
  rsv <- confidence_interval(ess.rsv.samp[[i]], .95)
  mean.asv[i] <- asv[1]
  mean.rsv[i] <- rsv[1]
  lower.asv[i] <- asv[2]
  lower.rsv[i] <- rsv[2]
  upper.asv[i] <- asv[3]
  upper.rsv[i] <- rsv[3]
}

plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "coral", main = "Two parallel chains", xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)))
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "blue")
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "coral")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"),col=c("coral", "blue"), lty=1, cex=1.2)
dev.off()


#####################################################################
#####################################################################
########## log-Frobenius norm density plots and coverage prob #######
#####################################################################
#####################################################################

########################################
####### A=1, B=3, C=8 ##################
########################################

A <- 1
B <- 3
C <- 8
check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5)

##### two chains

m <- 2
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))
  
  pdf(file = paste(paste("Out/", A, B, C, "/density_m", m, "_n", nsim, sep = "_"), ".pdf", sep = ""))
  plot(density(log(apply(asv.samp, 3, norm, type = "F"))), col = "coral", main = "", xlab = "log-frobenius norm",
       xlim = range(log(apply(asv.samp, 3, norm, type = "F")), log(apply(rsv.samp, 3, norm, type = "F"))))
  lines(density(log(apply(rsv.samp, 3, norm, type = "F"))), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

##### five chains

m <- 5
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))
  
  pdf(file = paste(paste("Out/", A, B, C, "/density_m", m, "_n", nsim, sep = "_"), ".pdf", sep = ""))
  plot(density(log(apply(asv.samp, 3, norm, type = "F"))), col="coral", main = "", xlab = "log-frobenius norm",
       xlim = range(log(apply(asv.samp, 3, norm, type = "F")), log(apply(rsv.samp, 3, norm, type = "F"))))
  lines(density(log(apply(rsv.samp, 3, norm, type = "F"))), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}


########################################
####### A=1, B=10, C=7 #################
########################################

A <- 1
B <- 10
C <- 7
check.pts <- c(1e3, 2e3, 5e3, 1e4, 5e4)

##### two chains

m <- 2
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))
  
  pdf(file = paste(paste("Out/", A, B, C, "/density_m", m, "_n", nsim, sep = "_"), ".pdf", sep = ""))
  plot(density(log(apply(asv.samp, 3, norm, type = "F"))), col="coral", main = "", xlab = "log-frobenius norm",
       xlim = range(log(apply(asv.samp, 3, norm, type = "F")), log(apply(rsv.samp, 3, norm, type = "F"))))
  lines(density(log(apply(rsv.samp, 3, norm, type = "F"))), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

##### five chains

m <- 5
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))
  
  pdf(file = paste(paste("Out/", A, B, C, "/density_m", m, "_n", nsim, sep = "_"), ".pdf", sep = ""))
  plot(density(log(apply(asv.samp, 3, norm, type = "F"))), col="coral", main = "", xlab = "log-frobenius norm",
       xlim = range(log(apply(asv.samp, 3, norm, type = "F")), log(apply(rsv.samp, 3, norm, type = "F"))))
  lines(density(log(apply(rsv.samp, 3, norm, type = "F"))), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

