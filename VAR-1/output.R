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

pdf(file = "Out/var-frob.pdf", height = 5, width = 5)
plot(conv.pts,log(apply(asv.samp, 3, norm, type = "F")), type = "l", col = "darkorange", lwd = 2,
     main = "", xlab = "Simulation size", ylab = "Variance",
     ylim = range(c(log(apply(rsv.samp, 3, norm, type = "F"))), log(norm(truth, type = "F"))))
lines(conv.pts, log(apply(rsv.samp, 3, norm, type = "F")), col="royalblue", lwd = 2)
abline(h = log(norm(truth, type = "F")), col = "green3", lwd = 2)
legend("bottomright", legend=c("ASV", "RSV", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, lwd = 2)
dev.off()

###########################################################
######### Effective Sample Size running plot##############
##########################################################

load(file = paste(paste("Out/run_data", m, min, max, sep = "_"), ".Rdata", sep = ""))

ind <- seq(1, length(conv.pts), by = 2)

a <- lapply(ess.asv.samp, log)
r <- lapply(ess.rsv.samp, log)
se.a <- 2*Reduce("rbind", lapply(a, sd))/sqrt(length(a[[1]]))
se.r <- 2*Reduce("rbind", lapply(r, sd))/sqrt(length(r[[1]]))
a <- Reduce("rbind", lapply(a, mean))
r <- Reduce("rbind", lapply(r, mean))

pdf(file = "Out/var-ess.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts[ind], y0 = (a - se.a)[ind], y1 = (a + se.a)[ind], col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts[ind], y0 = (r - se.r)[ind], y1 = (r + se.r)[ind], col = adjustcolor("royalblue", alpha.f = .50))
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
     ylab = "Y component", main = "", col = "darkorange")
points(mc.chain.list[[5]][1:nsim,], col = "royalblue")
legend("bottomright", legend = c("Chain-1", "Chain-5"), col = c("darkorange", "royalblue"), pch = 19)
dev.off()


nsim <- 1e5


pdf(file = ("Out/var-sp_n1e4.pdf"), height = 5, width = 5)
plot(mc.chain.list[[1]][1:nsim,], xlim = c(-(m/2)*sqrt(diag(target)[1]), (m/2)*sqrt(diag(target)[1])),
     ylim = c(-(m/2)*sqrt(diag(target)[2]), (m/2)*sqrt(diag(target)[2])), xlab = "X component",
     ylab = "Y component", main = "", col = "darkorange")
points(mc.chain.list[[5]][1:nsim,], col = "royalblue")
legend("bottomright", legend = c("Chain-1", "Chain-5"), col = c("darkorange", "royalblue"), pch = 19)
dev.off()

######################################################
###################### coverage prob #################
######################################################

for (j in 1:length(check.pts)){

  nsim <- check.pts[j]
  load(file = paste(paste("Out/out", m, nsim, sep = "_"), ".Rdata", sep = ""))

  print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, " are: "))
  print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))

}
