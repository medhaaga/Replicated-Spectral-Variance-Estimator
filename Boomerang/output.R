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


########################################################
#####Visualization through perspective plots############
########################################################

#Perspective plots for both simulation settings

A <- 1
B <- 3
C <- 8

pdf(file = paste("Out/", A, B, C, "/3d_density_plot.pdf", sep = "_"), height = 5, width = 5)
samples <- perspective(A, B, C, 100)
persp(samples$x, samples$y, samples$z, main="", col = "white", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()

A <- 1
B <- 10
C <- 7

pdf(file = paste("Out/", A, B, C, "/3d_density_plot.pdf", sep = "_"), height = 5, width = 5)
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

chain <- array(0, dim = c(1e3, p, m))
chain[,,1] <- markov.chain(A, B, C, 1e3, start[1,])
chain[,,2] <- markov.chain(A, B, C, 1e3, start[2,])

pdf(file = paste(paste("Out/", A, B, C, "/boom-sp", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(chain[,,1], col = "darkorange", xlim = range(chain[,1,]), ylim = range(chain[,2,]), xlab = "X-component", ylab = "Y-component", main = "")
points(chain[,,2], col = "royalblue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("darkorange", "royalblue"), pch=1)
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

chain <- array(0, dim = c(1e3, p, m))
chain[,,1] <- markov.chain(A, B, C, 1e3, start[1,])
chain[,,2] <- markov.chain(A, B, C, 1e3, start[2,])

pdf(file = paste(paste("Out/", A, B, C, "/boom-sp", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(chain[,,1], col = "darkorange", xlim = range(chain[,1,]), ylim = range(chain[,2,]),
     xlab = "X-component", ylab = "Y-component", main = "")
points(chain[,,2], col = "royalblue")
legend("topright", legend=c("chain-1", "chain-2"),col=c("darkorange", "royalblue"), pch = 1)
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
##### five chains ######

m <- 5
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/boom-frob", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "darkorange", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "log Frobenius norm",
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="royalblue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
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
##### five chains ######

m <- 5
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("Out/", A, B, C, "/boom-frob", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "darkorange", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "log Frobenius norm",
     ylim = range(log(apply(rsv.samp, 3, norm, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="royalblue", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
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
conv.pts <- seq(min, max, step)

#### five chains

m <- 5
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

a <- lapply(ess.asv.samp, log)
r <- lapply(ess.rsv.samp, log)
se.a <- 2*Reduce("rbind", lapply(a, sd))/sqrt(length(a[[1]]))
se.r <- 2*Reduce("rbind", lapply(r, sd))/sqrt(length(r[[1]]))
a <- Reduce("rbind", lapply(a, mean))
r <- Reduce("rbind", lapply(r, mean))

pdf(file = paste(paste("Out/", A, B, C, "/boom-ess", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1, lwd = 2)
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


#### five chains

m <- 5
load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

a <- lapply(ess.asv.samp, log)
r <- lapply(ess.rsv.samp, log)
se.a <- 2*Reduce("rbind", lapply(a, sd))/sqrt(length(a[[1]]))
se.r <- 2*Reduce("rbind", lapply(r, sd))/sqrt(length(r[[1]]))
a <- Reduce("rbind", lapply(a, mean))
r <- Reduce("rbind", lapply(r, mean))

pdf(file = paste(paste("Out/", A, B, C, "/boom-ess", A, B, C, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
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
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

##### five chains

m <- 5
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))

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
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

##### five chains

m <- 5
for (i in 1:length(check.pts)){
  nsim <- check.pts[i]
  load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"),".Rdata", sep = ""))
  print(paste("Coverage probabilities for nsim =  ", nsim, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
}

