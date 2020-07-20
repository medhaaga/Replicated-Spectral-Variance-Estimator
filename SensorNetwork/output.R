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

###################################################
### Data 
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

#######################################
####Model Parameters#################
######################################

m <- 5
p <- 8
j.scale <- rep(0.5, 4)

####### Starting Values####################

start1 <- c(-0.1, 0.5, -0.1, -0.2, 0.1, 0.1, -0.5, -0.5)
aux1 <- runif(n=8, min=min(start1), max=max(start1))
start2 <- c(0.0, 0.6, 0.1, 0.1, 0.2, 0.2, 1.0, 0.0)
aux2 <- runif(n=8, min=min(start2), max=max(start2))
start3 <- c(0.2, 0.7, 0.5, 0.4, 0.5, 0.3, 0.5, 0.5)
aux3 <- runif(n=8, min=min(start3), max=max(start3))
start4 <- c(0.4, 0.8, 0.8, 0.6, 0.7, 0.4, 1.0, 1.0)
aux4 <- runif(n=8, min=min(start4), max=max(start4))
start5 <- c(0.7, 1.0, 1.2, 0.9, 0.9, 0.5, 1.5, 1.5)
aux5 <- runif(n=8, min=min(start5), max=max(start5))

start <- rbind(start1, start2, start3, start4, start5)
aux <- rbind(aux1, aux2, aux3, aux4, aux5)

min <- 500
max <- 2e5
step <- 500
conv.pts <- seq(min, max, step)

####################################################
####################################################
###### log-frobenius norm running plots#############
#####################################################
#####################################################

load(file = paste("Out/conv_data_m", m, "_min", min, "_max", max, ".Rdata", sep = ""))

mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))


for (i in 1:length(conv.pts)){
  vec.asv <- rep(0, rep)
  vec.rsv <- rep(0, rep)
  for (j in 1:rep){
    vec.asv[j] <- log(norm(asv[[j]][,,i], type = "F"))
    vec.rsv[j] <- log(norm(rsv[[j]][,,i], type = "F"))
  }
  asv.ci <- confidence_interval(vec.asv, .95)
  rsv.ci <- confidence_interval(vec.rsv, .95)
  mean.asv[i] <- asv.ci[1]
  mean.rsv[i] <- rsv.ci[1]
  lower.asv[i] <- asv.ci[2]
  lower.rsv[i] <- rsv.ci[2]
  upper.asv[i] <- asv.ci[3]
  upper.rsv[i] <- rsv.ci[3]
}


pdf(file = paste("Out/run_plot-Frob.pdf", sep = ""))
plot(conv.pts, mean.asv, type = "l", col = "orange", main = "", lwd = 2,
     xlab = "Simulation size", ylab = "log-Frobenius norm", ylim = range(mean.asv, mean.rsv))
lines(conv.pts, mean.rsv, col="navy", lwd = 2)
segments(x0 = conv.pts, y0 = lower.asv, y1 = upper.asv, col = "orange")
segments(x0 = conv.pts, y0 = lower.rsv, y1 = upper.rsv, col = "navy")
legend("bottomright", legend=c("ASV", "RSV"),col=c("orange", "navy"), lty=1, lwd=2, cex=1)
dev.off()

##############################################
##############################################
####### Effective Sample Size ################
##############################################
##############################################


mean.asv <- rep(0, length(conv.pts))
mean.rsv <- rep(0, length(conv.pts))
lower.asv <- rep(0, length(conv.pts))
lower.rsv <- rep(0, length(conv.pts))
upper.asv <- rep(0, length(conv.pts))
upper.rsv <- rep(0, length(conv.pts))


for (i in 1:length(conv.pts)){
  vec.asv <- rep(0, rep)
  vec.rsv <- rep(0, rep)
  for (j in 1:rep){
    vec.asv[j] <- log(ess.asv[[j]][i])
    vec.rsv[j] <- log(ess.rsv[[j]][i])
  }
  asv.ci <- confidence_interval(vec.asv, .95)
  rsv.ci <- confidence_interval(vec.rsv, .95)
  mean.asv[i] <- asv.ci[1]
  mean.rsv[i] <- rsv.ci[1]
  lower.asv[i] <- asv.ci[2]
  lower.rsv[i] <- rsv.ci[2]
  upper.asv[i] <- asv.ci[3]
  upper.rsv[i] <- rsv.ci[3]
}


pdf(file = paste("Out/run_plot-ess.pdf", sep = ""))
plot(conv.pts, mean.asv, type = "l", col = "orange", main = "", lwd = 2,
     xlab = "Simulation size", ylab = "log-Frobenius norm", ylim = range(mean.asv, mean.rsv))
lines(conv.pts, mean.rsv, col="navy", lwd = 2)
segments(x0 = conv.pts, y0 = lower.asv, y1 = upper.asv, col = "orange")
segments(x0 = conv.pts, y0 = lower.rsv, y1 = upper.rsv, col = "navy")
legend("topright", legend=c("ASV", "RSV"),col=c("orange", "navy"), lty=1, lwd=2, cex=1)
dev.off()

################################################
################################################
### Scatter plots for four unknown locations ###
################################################
################################################

load(file = "Out/two_chains.Rdata")

nsim1 <- 1e4
nsim2 <- 1e5

pdf(file = "Out/sp-loc1.pdf", height = 4)
par(mfrow = c(1, 2))

# X_1 for nsim1
plot(res.ram1$x[1:nsim1, c(1, 2)], pch = 1, xlim = range(res.ram1$x[1:nsim1,1], res.ram2$x[1:nsim1,1], 0.5748),
     ylim = range(res.ram1$x[1:nsim1,2], res.ram2$x[1:nsim1,2]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim1,c(1,2)], pch = 1, xlim = range(res.ram1$x[1:nsim1,1], res.ram2$x[1:nsim1,1]),
       ylim = range(res.ram1$x[1:nsim1,2], res.ram2$x[1:nsim1,2]), xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)

#X_1 for nsim2

plot(res.ram1$x[1:nsim2, c(1, 2)], pch = 1, xlim = range(res.ram1$x[1:nsim2,1], res.ram2$x[1:nsim2,1]),
     ylim = range(res.ram1$x[1:nsim2,2], res.ram2$x[1:nsim2,2]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim2,c(1,2)], pch = 1, xlim = range(res.ram1$x[1:nsim2,1], res.ram2$x[1:nsim2,1]),
       ylim = range(res.ram1$x[1:nsim2,2], res.ram2$x[1:nsim2,2]), xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)
dev.off()

pdf(file = "Out/sp-loc2.pdf", height = 4)

par(mfrow = c(1,2))
#X_2 for nsim1
plot(res.ram1$x[1:nsim1, c(3, 4)], pch = 1, xlim = range(res.ram1$x[1:nsim1,3], res.ram2$x[1:nsim1,3]),
     ylim = range(res.ram1$x[1:nsim1,4], res.ram2$x[1:nsim1,4]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim1,c(3,4)], pch = 1, xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)

#X_2 for nsim2

plot(res.ram1$x[1:nsim2, c(3, 4)], pch = 1, xlim = range(res.ram1$x[1:nsim2,3], res.ram2$x[1:nsim2,3]),
     ylim = range(res.ram1$x[1:nsim2,4], res.ram2$x[1:nsim2,4]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim2,c(3,4)], pch = 1, xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)

dev.off()

pdf(file = "Out/sp-loc3.pdf", height = 4)
par(mfrow = c(1,2))
#X_3 for nsim1

plot(res.ram1$x[1:nsim1, c(5, 6)], pch = 1, xlim = range(res.ram1$x[1:nsim1,5], res.ram2$x[1:nsim1,5]),
     ylim = range(res.ram1$x[1:nsim1,6], res.ram2$x[1:nsim1,6]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim1,c(5,6)], pch = 1, xlim = range(res.ram1$x[1:nsim1,5], res.ram2$x[1:nsim1,5]),
       ylim = range(res.ram1$x[1:nsim1,6], res.ram2$x[1:nsim1,6]), xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)


#X_3 for nsim2

plot(res.ram1$x[1:nsim2, c(5, 6)], pch = 1, xlim = range(res.ram1$x[1:nsim2,5], res.ram2$x[1:nsim2,5]),
     ylim = range(res.ram1$x[1:nsim2,6], res.ram2$x[1:nsim2,6]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim2,c(5,6)], pch = 1, xlim = range(res.ram1$x[1:nsim2,5], res.ram2$x[1:nsim2,5]),
       ylim = range(res.ram1$x[1:nsim2,6], res.ram2$x[1:nsim2,6]), xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)
dev.off()

pdf("Out/sp-loc4.pdf", height = 4)

par(mfrow = c(1,2))
#X_4 for nsim1

plot(res.ram1$x[1:nsim1, c(7, 8)], pch = 1, xlim = range(res.ram1$x[1:nsim1,7], res.ram2$x[1:nsim1,7]),
     ylim = range(res.ram1$x[1:nsim1,8], res.ram2$x[1:nsim1,8]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim1,c(7,8)], pch = 1, xlim = range(res.ram1$x[1:nsim1,7], res.ram2$x[1:nsim1,7]),
       ylim = range(res.ram1$x[1:nsim1,8], res.ram2$x[1:nsim1,8]), xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)

#X_4 for nsim2

plot(res.ram1$x[1:nsim2, c(7, 8)], pch = 1, xlim = range(res.ram1$x[1:nsim2,7], res.ram2$x[1:nsim2,7]),
     ylim = range(res.ram1$x[1:nsim2,8], res.ram2$x[1:nsim2,8]), xlab = "", ylab = "", main = "", col = "orange")
points(res.ram2$x[1:nsim2,c(7,8)], pch = 1, xlab = "", ylab = "", main = "", col = "navy")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)
dev.off()

########################################################
########################################################
### Trace for x-coordinate of four unknown locations ###
########################################################
########################################################

load(file = "Out/two_chains.Rdata")

pdf(file = "Out/trace.pdf")
par(mfrow = c(2,2))

plot.ts(res.ram1$x[,1])
plot.ts(res.ram1$x[,3])
plot.ts(res.ram1$x[,5])
plot.ts(res.ram1$x[,7])
dev.off()