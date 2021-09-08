set.seed(1)
library(fields)
library(graphics)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
sourceCpp("lag.cpp")
source("functions.R")


m <- 5
p <- 100
rep <- 10
min <- 5e2
max <- 5e4
step <- 500
conv.pts <- seq(min, max, step)


#####################################################
##### Setting 1: Good mixing convergence plots ######
#####################################################

set <- 1

############## log Frobenius norm running plots ###########

load( file = "Out/var-set1_truth.Rdata")
load(file = paste("Out/conv_data_min", min, "_max", max, "_set", set, ".Rdata", sep = ""))

a <- lapply(asv, function(x) x/norm(truth, "F"))
r <- lapply(rsv, function(x) x/norm(truth, "F") )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = "Out/var-set1_frob.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a, 1)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = 1, col = "green3", lwd = 2)
legend("bottomright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, lwd = 2)
dev.off()


######## Effective Sample Size running plot##############

load(file = paste("Out/conv_data_min", min, "_max", max, "_set", set, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = "Out/var-set1_ess.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r, log((det(target)/det(truth))^(1/p))), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = log((det(target)/det(truth))^(1/p)), col = "green3", lwd = 2)
legend("topright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, cex=1, lwd = 2)
dev.off()

#####################################################
##### Setting 2: Bad mixing convergence plots ######
#####################################################

set <- 2

############## log Frobenius norm running plots ###########

load( file = "Out/var-set2_truth.Rdata")
load(file = paste("Out/conv_data_min", min, "_max", max, "_set", set, ".Rdata", sep = ""))

a <- lapply(asv, function(x) x/norm(truth, "F"))
r <- lapply(rsv, function(x) x/norm(truth, "F") )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = "Out/var-set2_frob.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a, 1)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = 1, col = "green3", lwd = 2)
legend("bottomright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, lwd = 2)
dev.off()

######### Effective Sample Size running plot##############

load(file = paste("Out/conv_data_min", min, "_max", max, "_set", set, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = "Out/var-set2_ess.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r, (1/p)*(sum(log(eigen(target)$values)) - sum(log(eigen(truth)$values)))), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = (1/p)*(sum(log(eigen(target)$values)) - sum(log(eigen(truth)$values))), col = "green3", lwd = 2)
legend("topright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, cex=1, lwd = 2)
dev.off()

