set.seed(1)
library(fields)library(cubature)
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


###############################################
##############Model parameters#################
###############################################

c0 <- 13
d0 <- 1
bars <- 6
p <- bars+1
m <- 2
min <- 5e2
max <- 1e5
step <- 500
conv.pts <- seq(min, max, step)

##############################################
##############################################
## log Frobenius norm Running plots ##########
##############################################
##############################################

load(file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

frob.asv <- log(apply(asv.samp, 3, norm, type = "F"))
frob.rsv <- log(apply(rsv.samp, 3, norm, type = "F"))
pdf(file = paste(paste("Out/run_plot-Frob", sep = "_"), ".pdf", sep = ""))
plot(conv.pts, frob.asv, type = "l", col="coral", main = "", xlab = "Simulation size", 
     ylab = "log-Frobenius norm", ylim = range(frob.asv, frob.rsv))
lines(conv.pts, frob.rsv, col="royalblue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=.75)
dev.off()

##################################################
##################################################
######### Effective Sample Size ##################
##################################################
##################################################

pdf(file = paste("Out/run_plot-ess.pdf", sep = "_"), height = 4)

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
plot(conv.pts[1:100], mean.asv[1:100], type = "l", col = "coral", main = "", 
     xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)), lwd=2)
lines(conv.pts[1:100], mean.rsv[1:100], type = "l", col = "royalblue", lwd=2)
segments(x0 = conv.pts[1:100], y0 = lower.asv[1:100], y1 = upper.asv[1:100], col = "red")
segments(x0 = conv.pts[1:100], y0 = lower.rsv[1:100], y1 = upper.rsv[1:100], col = "blue")
legend("topright", legend=c("ASV", "RSV"), col=c("coral", "royalblue"), lty=1, cex=1.2, lwd=2)
dev.off()

########################################
########################################
###########trace plots##################
########################################
########################################


chain <- matrix(MCMCpoissonChange(mida~ 1, m = bars, c0=c0, d0=d0, 
                                       marginal.likelihood="none", mcmc = 1e4, burnin = 0, thin = 1, seed = 500), nrow = 1e4, ncol = p)

for(i in 1:p){
  pdf(file = paste(paste("Out/trace_component", i, sep = "_"), ".pdf", sep = ""), height = 4)
  plot.ts(chain[,i], xlab = "Time", ylab = paste("Component-", j, sep = ""), main = "")
  dev.off()
}

