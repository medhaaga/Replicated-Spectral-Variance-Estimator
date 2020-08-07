set.seed(1)
library(cubature)
library(graphics)
library(pracma)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
library(plotly)
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

check.pts <- c(5e3, 1e4, 5e4, 1e5)
r <- length(check.pts)
freq <- 1e2
c.prob <- .95
min <- 500
max <- 2e5
step <- 500
conv.pts <- seq(min, max, step)


p <- 8
load(file = paste("Out/conv_data_m", m, "_min", min, "_max", max, ".Rdata", sep = ""))

########################################
########################################
### log Frobenius norm running plot ####
########################################
########################################

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)
pdf(file = paste("Out/sensor-frob.pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts[ind], y0 = (a - se.a)[ind], y1 = (a + se.a)[ind], col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts[ind], y0 = (r - se.r)[ind], y1 = (r + se.r)[ind], col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

########################################
########################################
########### log ESS running plot #######
########################################
########################################

pdf(file = paste("Out/sensor-ess.pdf", sep = "_"), height = 5, width = 5)

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r, log(minESS(p, alpha = .1, eps = .1)/(m*conv.pts))), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
lines(conv.pts, log(minESS(p, alpha = .1, eps = 0.1)/(m*conv.pts)), col = "green3", lwd=2)
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1.2, lwd=2)
dev.off()

###############################################
###############################################
# Scatter plots
###############################################
###############################################


load(file = "Out/five_chains.Rdata")
nsim1 = 1e5

### Scatter plot - Location-1

df <- as.data.frame(mc.chain.list[[1]][,c(1,2)])
names(df) <- c("X", "Y")
fig <- ggplot(df, aes(X, Y)) + geom_density_2d_filled() + theme_minimal() + scale_fill_brewer(palette="Blues") + theme(legend.position="none") + xlim(-0.25,0.75)
pdf("Out/sensor-contour_loc1.pdf", height = 5, width = 5)
print(fig)
dev.off()

### Trace plot - Location-1

pdf(file = "Out/sensor-trace_loc1.pdf", height = 5, width = 5)
plot.ts(mc.chain.list[[1]][1:nsim1,1], ylim = range(mc.chain.list[[1]][,1], mc.chain.list[[5]][,1]), ylab = expression(x[11]), col = rgb(8, 69, 148, maxColorValue = 255))
par(new = TRUE)
plot.ts(mc.chain.list[[5]][1:nsim1,1], ylim = range(mc.chain.list[[1]][,1], mc.chain.list[[5]][,1]), yaxt='n', xaxt='n', ylab = "", col = rgb(107, 174, 214, maxColorValue = 255))
dev.off()
