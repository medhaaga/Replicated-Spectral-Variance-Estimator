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
sourceCpp("lag.cpp")
source("functions.R")

#################################################
##### Model parameters for both setttings #######
#################################################

p <- 2

A1 <- 1
B1 <- 3
C1 <- 8

A2 <- 1
B2 <- 10
C2 <- 7

########################################################
#####Visualization through perspective plots############
########################################################

#### Setting-1

pdf(file = paste(paste("Out/", A1, B1, C1, "/3d_density_plot", A1, B1, C1, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
samples <- perspective(A1, B1, C1, 100)
persp(samples$x, samples$y, samples$z, main="", col = "white", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()

#### Setting-2

pdf(file = paste(paste("Out/", A2, B2, C2, "/3d_density_plot", A2, B2, C2, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
samples <- perspective(A2, B2, C2, 200)
persp(samples$x, samples$y, samples$z, main="", col = "white", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()


########################################################
########################################################
####### Scatter plots for both simulation settings######
########################################################
########################################################

########### Setting-1 ##############

set.seed(1)  # use seed for reproducible results
m <- 2

start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C1*(2^(2-i)))
  start[m-i+1,] <- c(C1*(2^(2-i)), 0)
}

chain <- array(0, dim = c(1e3, p, m))
chain[,,1] <- markov.chain(A1, B1, C1, 1e3, start[1,])
chain[,,2] <- markov.chain(A1, B1, C1, 1e3, start[2,])

pdf(file = paste(paste("Out/", A1, B1, C1, "/boom-sp", A1, B1, C1, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(chain[,,1], col = 'dodgerblue4', xlim = range(chain[,1,]), ylim = range(chain[,2,]), xlab = "X-component", ylab = "Y-component", main = "")
points(chain[,,2], col = "steelblue1")
legend("topright", legend=c("chain-1", "chain-2"),col=c("dodgerblue4", "steelblue1"), pch=1, cex=1.2)
dev.off()


########### Setting-2 ##############

set.seed(1)  # use seed for reproducible results
m <- 2

start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C2*(2^(2-i)))
  start[m-i+1,] <- c(C2*(2^(2-i)), 0)
}

chain <- array(0, dim = c(1e3, p, m))
chain[,,1] <- markov.chain(A2, B2, C2, 1e3, start[1,])
chain[,,2] <- markov.chain(A2, B2, C2, 1e3, start[2,])

pdf(file = paste(paste("Out/", A2, B2, C2, "/boom-sp", A2, B2, C2, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(chain[,,1], col = 'dodgerblue4', xlim = range(chain[,1,]), ylim = range(chain[,2,]), xlab = "X-component", ylab = "Y-component", main = "")
points(chain[,,2], col = "steelblue1")
legend("topright", legend=c("chain-1", "chain-2"),col=c("dodgerblue4", "steelblue1"), pch=1, cex=1.2)
dev.off()


###################################################################
###################################################################
########## Running plots for log Frobenius norm ###################
###################################################################
###################################################################


##########################################
######### Setting-1 ##################
##########################################


min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)

##### five chains ######

m <- 5
load(file = paste(paste("Out/", A1, B1, C1, "/conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = paste(paste("Out/", A1, B1, C1, "/boom-frob", A1, B1, C1, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

###########################################
######### Setting-2 ##################
###########################################

min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)

##### five chains ######

m <- 5
load(file = paste(paste("Out/", A2, B2, C2, "/conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = paste(paste("Out/", A2, B2, C2, "/boom-frob", A2, B2, C2, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

##########################################################
##########################################################
#######Running plots for Effective Sample Size############
##########################################################
##########################################################


###########################################
############# Setting-1 ###############
###########################################


min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)

#### five chains

m <- 5
load(file = paste(paste("Out/", A1, B1, C1, "/conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = paste(paste("Out/", A1, B1, C1, "/boom-ess", A1, B1, C1, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1, lwd = 2)
dev.off()


#############################################
############# Setting-2 ###############
############################################


min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)

#### five chains

m <- 5
load(file = paste(paste("Out/", A2, B2, C2, "/conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = paste(paste("Out/", A2, B2, C2, "/boom-ess", A2, B2, C2, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1, lwd = 2)
dev.off()

#########################################
#########################################
##########  coverage prob ###############
#########################################
#########################################

########################################
####### Setting-1##################
########################################

check.pts <- c(1e3, 5e3, 1e4, 2e4, 5e4, 1e5)
freq <- 1000
##### two chains

m <- 2
load(file = paste(paste("Out/", A1, B1, C1, "/out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  
}


##### five chains

m <- 5
load(file = paste(paste("Out/", A1, B1, C1, "/out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  
}


########################################
####### Setting-2 #################
########################################

check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4, 5e4)
freq <- 1000

##### two chains

m <- 2

load(file = paste(paste("Out/", A2, B2, C2, "/out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  
}


##### five chains

m <- 5
load(file = paste(paste("Out/", A2, B2, C2, "/out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  
}
