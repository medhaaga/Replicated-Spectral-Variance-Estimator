set.seed(1)
library(fields)
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
library(ergm)


#############Data Modification#####################
data(faux.magnolia.high)
magnolia = faux.magnolia.high

#Ensure well connected component only 
notwellconnected = which(component.largest(magnolia, connected="weak")==FALSE)
delete.vertices(magnolia, notwellconnected)
dp = dataPrepHighSchool(magnolia)

###################################################
### Model parameters

start <- c(279, 360)
p <- 5
m <- 2   #no of parallel chains
min <- 1e2
max <- 1e4
step <- 100
conv.pts <- seq(min, max, step)

###################################################
###################################################
######log Frobenius norm running plots ###########
###################################################
###################################################

load(file = paste(paste("Out/conv_data", min, max, sep = "_"), ".Rdata", sep = ""))

frob.asv <- log(apply(asv.samp[3:5, 3:5,], 3, norm, type = "F"))
frob.rsv <- log(apply(rsv.samp[3:5, 3:5,], 3, norm, type = "F"))

pdf(file = paste(paste("Out/run_plot-Frob", sep = "_"), ".pdf", sep = ""))
plot(conv.pts, frob.asv, type = "l", col="coral", main = "", xlab = "Simulation size", 
     ylab = "log-Frobenius norm", ylim = range(frob.asv, frob.rsv), lwd=2)
lines(conv.pts, frob.rsv, col="royalblue", lwd=2)
legend("topright", legend = c("ASV", "RSV"), col=c("red", "blue"), lty=1, cex=1, lwd=2)
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

plot(conv.pts, mean.asv, type = "l", col = "coral", main = "", 
     xlab = "Simulation size", ylab = "ESS/mn", ylim = range(c(mean.asv, mean.rsv)), lwd=2)
lines(conv.pts, mean.rsv, type = "l", col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = lower.asv, y1 = upper.asv, col = "red")
segments(x0 = conv.pts, y0 = lower.rsv, y1 = upper.rsv, col = "blue")
legend("topright", legend=c("ASV", "RSV"), col=c("coral", "royalblue"), lty=1, cex=1.2, lwd=2)
dev.off()

