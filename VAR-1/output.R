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

phi <- diag(c(.999, .001))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

############3Simuation settings###########

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
check.pts <- c(5e2, 1e3, 5e3, 1e4, 5e4, 1e5)
rep <- 50
freq <- 1000
min <- 5e2
max <- 5e4
step <- 500
conv.pts <- seq(min, max, step)

###########################################################
############## log Frobenius norm running plots ###########
###########################################################

load(file = paste("Out/conv_data_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = "Out/var-frob.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a, log(norm(truth, type = "F")))), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = log(norm(truth, type = "F")), col = "green3", lwd = 2)
legend("bottomright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, lwd = 2)
dev.off()

###########################################################
######### Effective Sample Size running plot##############
##########################################################

load(file = paste("Out/conv_data_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = "Out/var-ess.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = log((det(target)/det(truth))^(1/p)), col = "green3", lwd = 2)
legend("topright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, cex=1, lwd = 2)
dev.off()


########################################3
#########################################
####  Scatter plots #####################
#########################################
#########################################

load(file = "Out/five_chains.Rdata")

nsim <- 1e3

pdf(file = ("Out/var-sp_n1e3.pdf"), height = 5, width = 5)
plot(mc.chain.list[[1]][1:nsim,], xlim = c(-(m/2)*sqrt(diag(target)[1]), (m/2)*sqrt(diag(target)[1])),
     ylim = c(-(m/2)*sqrt(diag(target)[2]), (m/2)*sqrt(diag(target)[2])), xlab = "X component",
     ylab = "Y component", main = "", col = "dodgerblue4")
points(mc.chain.list[[5]][1:nsim,], col = "steelblue1")
legend("bottomright", legend = c("Chain-1", "Chain-5"), col = c("dodgerblue4", "steelblue1"), pch = 19)
dev.off()


nsim <- 1e5


pdf(file = ("Out/var-sp_n1e4.pdf"), height = 5, width = 5)
plot(mc.chain.list[[1]][1:nsim,], xlim = c(-(m/2)*sqrt(diag(target)[1]), (m/2)*sqrt(diag(target)[1])),
     ylim = c(-(m/2)*sqrt(diag(target)[2]), (m/2)*sqrt(diag(target)[2])), xlab = "X component",
     ylab = "Y component", main = "", col = "dodgerblue4")
points(mc.chain.list[[5]][1:nsim,], col = "steelblue1")
legend("bottomright", legend = c("Chain-1", "Chain-5"), col = c("dodgerblue4", "steelblue1"), pch = 19)
dev.off()

######################################################
###################### coverage prob #################
######################################################

load(file = paste("Out/out_check.pts_freq", freq, ".Rdata", sep = ""))

for (j in 1:length(check.pts)){

  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))

}
