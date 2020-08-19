library(rep.acf.ccf)

##########################################
##########################################
####### Ex1: Gaussian Mixtures ###########
##########################################
##########################################

log.density <- function(x, p, mu1, mu2, sd1, sd2){
  return(log(p*dnorm(x, mean = mu1, sd = sd1) + (1-p)*dnorm(x, mean = mu2, sd= sd2)))
}

m <- 2
p <- 0.7
mu1 <- -5
mu2 <- 5
sd1 <- 1
sd2 <- 0.5

load(file = "AllOut/gaussian-two_chains.Rdata")
chain1 <- mc.chain.list[[1]]
chain2 <- mc.chain.list[[2]]

################### Figure 1 #######################

## Figure 1a

x <- seq(-10, 10, length = 1e3)
pdf(file = "AllOut/gaussian-TargetTrace_n1e4.pdf", height = 5, width = 6)
plot(x, 20000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x", ylab = "", ylim = c(-10000,6000), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1[1:1e4], y = seq(-1, -1e4, -1), col = "lightskyblue", xlab = "", ylab = "", type = "l", ylim = c(-1e4, 6e3), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2[1:1e4], y = seq(-1, -1e4, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.7, lwd=2)
dev.off()

#Figure 1b
pdf(file = "AllOut/gaussian-TargetTrace_n1e5.pdf", height = 5, width = 6)
plot(x, 200000*exp(log.density(x, p, mu1, mu2, sd1, sd2)), type = "l", lwd=2, xlab = "x", ylab = "", ylim = c(-1e5,6e4), xlim = range(chain1,chain2), yaxt = 'n')
par(new = TRUE)
plot(x = chain1, y = seq(-1, -1e5, -1), col = "lightskyblue", xlab = "", ylab = "", type = "l", ylim = c(-1e5, 6e4), , xlim = range(chain1,chain2), yaxt = 'n')
lines(x = chain2, y = seq(-1, -1e5, -1), col = "plum3")
mtext(side = 2, text = "Time", line = 1)
legend("topleft", legend=c("Target", "Chain-1", "Chain-2"),col=c("black", "lightskyblue", "plum3"), lty=1, cex=0.7, lwd=2)
dev.off()

#################### Figure 2 #######################

lag.max <- 50
nsim1 <- 1e4
nsim2 <- 5e4

x <- list()
y <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]], ncol=1)[1:nsim1,]
  y[[i]] <- as.matrix(mc.chain.list[[i]], ncol=1)[1:nsim2,]
  x[[i]] <- as.matrix(x[[i]])
  y[[i]] <- as.matrix(y[[i]])
}

local.acf1 <- acf(as.matrix(x[[1]]), type = "correlation", lag.max = lag.max, plot = FALSE)
local.acf2 <- acf(as.matrix(y[[1]]), type = "correlation", lag.max = lag.max, plot = FALSE)
global.acf1 <- globalACF(x, chains = c(1), component = 1, lag.max = lag.max, type = "correlation", avg = FALSE, graph = FALSE)[[1]]
global.acf2 <- globalACF(y, chains = c(1), component = 1, lag.max = lag.max, type = "correlation", avg = FALSE, graph = FALSE)[[1]]

#### Figure 2a
pdf(file = "AllOut/gaussian-acf_hist.pdf", width = 10, height= 4)
par(mfrow = c(1,2))
plot(as.matrix(local.acf1$acf), type = 'h', ylab = "Autocorrelation", xlab = "Lag")
lines(as.matrix(local.acf2$acf), type = 'l', col = "steelblue1", lwd = 2)
plot(as.matrix(global.acf1$acf), type = 'h', ylim = c(min(local.acf1$acf), 1), ylab = "Autocorrelation", xlab = "Lag")
lines(as.matrix(global.acf2$acf), type = 'l', col = "steelblue1", lwd = 2)
dev.off()

### Figure 2b
pdf(file = "AllOut/gaussian-acf_1e4.pdf", width = 10, height= 4)
par(mfrow = c(1,2))
l <- globalACF(x, chains = 0, component = 1, lag.max = lag.max, mean = "local", type = "correlation", col = "royalblue", leg = FALSE)
g <- globalACF(x, chains = 0, component = 1, lag.max = lag.max, mean = "global", type = "correlation", col = "darkorange", leg = FALSE)
dev.off()


##########################################
##########################################
#### Ex-2 VAR (1) ########################
##########################################
##########################################

##################Figure 3 ###############

load(file = "AllOut/var-five_chains.Rdata")

m <- 5
p <- 2
component <- 1
lag.max <- 40

### Figure 3a
ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, lag.max = lag.max, chains = c(2), graph = FALSE, avg = FALSE)[[1]]
local.acf <- acf(x[[2]][, component], lag.max = lag.max, type = "correlation", plot = FALSE)

pdf(file = paste("AllOut/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()

### Figure 3b

ncrop <- 1e4
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, lag.max = lag.max, chains = c(2), graph = FALSE, avg = FALSE)[[1]]
local.acf <- acf(x[[2]][, component], lag.max = lag.max, type = "correlation", plot = FALSE)

pdf(file = paste("AllOut/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()

##################### Figure 4 ######################

min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)
load(file = "AllOut/var-truth.Rdata")

### Figure 4a

load(file = paste("AllOut/var-conv_data_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = "AllOut/var-frob.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a, log(norm(truth, type = "F")))), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = log(norm(truth, type = "F")), col = "green3", lwd = 2)
legend("bottomright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, lwd = 2)
dev.off()

### Figure 3b

load(file = paste("AllOut/var-conv_data_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = "AllOut/var-ess.pdf", height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
abline(h = log((det(target)/det(truth))^(1/p)), col = "green3", lwd = 2)
legend("topright", legend=c("A-SVE", "G-SVE", "Truth"),col=c("darkorange", "royalblue", "green3"), lty=1, cex=1, lwd = 2)
dev.off()

################# Table 1 ######################

load(file = paste("AllOut/var-out_check.pts_freq", freq, ".Rdata", sep = ""))
check.pts <- c(5e2, 1e3, 5e3, 1e4, 5e4, 1e5)
freq <- 1e3

for (j in 1:length(check.pts)){
  
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  
}

##############################################
##############################################
########## Ex3 Boomerang #####################
##############################################
##############################################

library(plot3D)
p <- 2
A1 <- 1
B1 <- 3
C1 <- 8

A2 <- 1
B2 <- 10
C2 <- 7

################## Figure 5 #############################

### Figure 5a

load(file = paste(paste("AllOut/boom-two_chains_sp", A1, B1, C1, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("AllOut/boom-2d_density_plot", A1, B1, C1, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
contour2D(x=samples$x, y=samples$y, z=samples$z, colkey = FALSE, xlim=c(0,10), ylim=c(0,10))
points(rbind(chain[,,1], chain[,,2]), cex=.2, col = c(rep("black",1e3), rep("darkorange",1e3)))
dev.off()

### Figure 5b

load(file = paste(paste("AllOut/boom-two_chains_sp", A2, B2, C2, sep = "_"), ".Rdata", sep = ""))

pdf(file = paste(paste("AllOut/boom-2d_density_plot", A2, B2, C2, sep = "_"), ".pdf", sep = ""), height = 5, width = 5)
contour2D(x=samples$x, y=samples$y, z=samples$z, colkey = FALSE, xlim=c(0,10), ylim=c(0,10))
points(rbind(chain[,,1], chain[,,2]), cex=.2, col = c(rep("black",1e3), rep("darkorange",1e3)))
dev.off()

########################## Figure  6 #####################

load(file = "AllOut/boom-five_chains_1_3_8.Rdata")
m <- 5
component <- 1
lag.max <- 50
nsim1 <- 1000
nsim2 <- 1e5

x <- list()
y <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:nsim1,]
  y[[i]] <- mc.chain.list[[i]][1:nsim2,]
}

global.acf1 <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$"avg_ACF"
global.acf2 <- globalACF(y, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$"avg_ACF"
local.acf1 <- globalACF(x, type = "correlation", lag.max = lag.max, mean = "local", component = 1, graph = FALSE)$"avg_ACF"
local.acf2 <- globalACF(y, type = "correlation", lag.max = lag.max, mean = "local", component = 1, graph = FALSE)$"avg_ACF"

pdf(file = paste(paste("AllOut/boom-acf", A1, B1, C1, sep = "_"), ".pdf", sep = ""), width = 10, height = 4)
par(mfrow = c(1,2))
plot(local.acf1, xlab = "Lag", ylab = "Autocorrelation", main = "")
lines(local.acf2$acf, col = "steelblue1", lwd=2)
plot(global.acf1, xlab = "Lag", ylab = "Autocorrelation", main = "")
lines(global.acf2$acf, col = "steelblue1", lwd=2)
dev.off()


###################### Figure 7 ##########################

load(file = "AllOut/boom-five_chains_1_10_7.Rdata")
nsim <- 1000
x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$"avg_ACF"
local.acf <- globalACF(x, type = "correlation", lag.max = lag.max, mean = "local", component = 1, graph = FALSE)$"avg_ACF"

pdf(file = paste(paste("AllOut/boom-acf", A2, B2, C2, "n", sep = "_"),  nsim, ".pdf", sep = ""), width = 10, height = 4)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = "")
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = "")
dev.off()

####################### Figure 8 ########################

min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)

m <- 5

### Figure 8a
load(file = paste(paste("AllOut/boom1-conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(a)
r <- Reduce("+", r)/length(r)

pdf(file = paste(paste("AllOut/boom-ess", A1, B1, C1, "m", sep = "_"), m, ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1, lwd = 2)
dev.off()

### Figure 8b

min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)
load(file = paste(paste("AllOut/boom2-conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(a)
r <- Reduce("+", r)/length(r)

pdf(file = paste(paste("AllOut/boom-ess", A2, B2, C2, "m", sep = "_"), m, ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd=2)
lines(conv.pts, r, col = "royalblue", lwd=2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1, lwd = 2)
dev.off()

####################### Figure 9 ####################

#### Figure 9a
min <- 5e2
max <- 1e5
step <- 5e2
conv.pts <- seq(min, max, step)
m <- 5
load(file = paste(paste("AllOut/boom1-conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv[1:10], function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv[1:10], function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(10)
se.r <- apply(r, 2, sd)/sqrt(10)

a <- colMeans(a)
r <- colMeans(r)

pdf(file = paste(paste("AllOut/boom-frob", A1, B1, C1, "m", sep = "_"), m, ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

#### Figure 9b

min <- 5e2
max <- 5e4
step <- 5e2
conv.pts <- seq(min, max, step)
m <- 5
load(file = paste(paste("AllOut/boom2-conv_data_m", sep = "_"), m, "_min", min, "_max", max, ".Rdata", sep = ""))

a <- lapply(asv[1:10], function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv[1:10], function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(10)
se.r <- apply(r, 2, sd)/sqrt(10)

a <- colMeans(a)
r <- colMeans(r)

pdf(file = paste(paste("AllOut/boom-frob", A2, B2, C2, "m", sep = "_"), m, ".pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

##################### Table 2 ########################

check.pts <- c(1e3, 5e3, 1e4, 2e4, 5e4, 1e5)
freq <- 1000

##### two chains

m <- 2
load(file = paste(paste("AllOut/boom1-out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
}

##### five chains

m <- 5
load(file = paste(paste("AllOut/boom1-out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
  }

##################### Table 3 ########################

check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4, 5e4)
freq <- 1e3

##### two chains

m <- 2
load(file = paste(paste("AllOut/boom2-out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
}

##### five chains

m <- 5
load(file = paste(paste("AllOut/boom2-out_check.pts_m", sep = "_"), m, "_freq", freq, ".Rdata", sep=""))

for (j in 1:length(check.pts)){
  nsim <- check.pts[j]
  print(paste("Coverage probabilities for nsim =  ", nsim))
  print(paste("ASV : ", mean(asv.coverage[[j]]), "GSV: ", mean(gsv.coverage[[j]])))
}

######################################################
######################################################
############## Ex3 Sensor Metwork ####################
######################################################
######################################################

p <- 8
m <- 5

######################## Figure 10 ########################

load(file = "AllOut/sensor-five_chains.Rdata")
nsim1 = 1e5

### Trace plot - Location-1

pdf(file = "AllOut/sensor-trace_loc1.pdf", height = 4, width = 10)
par(mfrow = c(1,2))
plot.ts(mc.chain.list[[1]][1:nsim1,1], ylim = range(mc.chain.list[[1]][,1], mc.chain.list[[5]][,1]), ylab = expression(x[11]), col = rgb(8, 69, 148, maxColorValue = 255))
par(new = TRUE)
plot.ts(mc.chain.list[[5]][1:nsim1,1], ylim = range(mc.chain.list[[1]][,1], mc.chain.list[[5]][,1]), yaxt='n', xaxt='n', ylab = "", col = rgb(107, 174, 214, maxColorValue = 255))

plot.ts(mc.chain.list[[1]][1:nsim1,2], ylim = range(mc.chain.list[[1]][,2], mc.chain.list[[5]][,2]), ylab = expression(x[12]), col = rgb(8, 69, 148, maxColorValue = 255))
par(new = TRUE)
plot.ts(mc.chain.list[[5]][1:nsim1,2], ylim = range(mc.chain.list[[1]][,2], mc.chain.list[[5]][,2]), yaxt='n', xaxt='n', ylab = "", col = rgb(107, 174, 214, maxColorValue = 255))
dev.off()

######################## Figure 11 ############################

load(file = "AllOut/sensor-five_chains.Rdata")
component <- 1
lag.max <- 40

### Figure 11a
ncrop <- 5e3

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

pdf(file = "AllOut/sensor-acf_n5e3.pdf", height = 4, width = 4)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
par(new = TRUE)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

### Figure 11b

ncrop <- 5e4

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

pdf(file = "AllOut/sensor-acf_n5e4.pdf", height = 4, width = 4)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
par(new = TRUE)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

######################## Figure 12 ###########################

min <- 500
max <- 2e5
step <- 500
conv.pts <- seq(min, max, step)

load(file = paste("AllOut/sensor-conv_data_m", m, "_min", min, "_max", max, ".Rdata", sep = ""))


### Figure 12a

a <- lapply(asv, function(x) log(apply(x, 3, norm, "F") ) )
r <- lapply(rsv, function(x) log(apply(x, 3, norm, "F") ) )
a <- Reduce("rbind", a)
r <- Reduce("rbind", r)

se.a <- apply(a, 2, sd)/sqrt(length(asv))
se.r <- apply(r, 2, sd)/sqrt(length(rsv))

a <- colMeans(a)
r <- colMeans(r)

pdf(file = paste("AllOut/sensor-frob.pdf", sep = ""), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "Log of Frobenium norm", ylim = range(c(a, r, r + se.r, a - se.a)), lwd = 2)
lines(conv.pts, r, col="royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("bottomright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, lwd = 2)
dev.off()

### Figure 12b

a <- lapply(ess.asv, log)
r <- lapply(ess.rsv, log)
se.a <- 2*apply(Reduce("rbind", a), 2, sd)/sqrt(length(a))
se.r <- 2*apply(Reduce("rbind", r), 2, sd)/sqrt(length(r))
a <- Reduce("+", a)/length(ess.asv)
r <- Reduce("+", r)/length(ess.rsv)

pdf(file = paste("AllOut/sensor-ess.pdf", sep = "_"), height = 5, width = 5)
plot(conv.pts, a, type = "l", col = "darkorange", main = "", xlab = "Simulation size", ylab = "log(ESS/mn)", ylim = range(a, r), lwd = 2)
lines(conv.pts, r, col = "royalblue", lwd = 2)
segments(x0 = conv.pts, y0 = (a - se.a), y1 = (a + se.a), col = adjustcolor("darkorange", alpha.f = .50))
segments(x0 = conv.pts, y0 = (r - se.r), y1 = (r + se.r), col = adjustcolor("royalblue", alpha.f = .50))
legend("topright", legend=c("A-SVE", "G-SVE"),col=c("darkorange", "royalblue"), lty=1, cex=1.2, lwd=2)
dev.off()


###############################################
###############################################
############ Ex4 Poisson Change ###############
###############################################
###############################################

################### Figure 13 #################

load(file = "AllOut/poisson-two_chains.Rdata")

n <- 5e4
pdf(file = paste("AllOut/poisson-trace_n", n, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot.ts(mc.chain.list[[1]][1:n,2], col = "steelblue1", xlab = "Time", ylab = "Component-2", main = "")
lines(1:n, mc.chain.list[[2]][1:n, 2], col = "dodgerblue4", xlab = "Time", ylab = "Component-2", main = "")
plot.ts(mc.chain.list[[1]][1:n,3], col = "steelblue1", xlab = "Time", ylab = "Component-3", main = "")
lines(1:n, mc.chain.list[[2]][1:n,3], col = "dodgerblue4", xlab = "Time", ylab = "Component-3", main = "")
dev.off()

################### Figure 14 ########################

load(file = "AllOut/poisson-two_chains.Rdata")
m <- 2
component <- 2

### Figure 14a

ncrop <- 1e3
x <- list()

for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("AllOut/poisson-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

### Figure 14b

ncrop <- 1e4
x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("AllOut/poisson-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

############################################
##########################################
######## Ex5 Network Crawling ##############
#############################################
#############################################

##################  Figure 15 #################

m <- 2
component <- 3
load(file = "AllOut/magnolia-two_chains.Rdata")

### Figure 15a

ncrop <- 1e2
x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("AllOut/magnolia-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

### Figure 15b

ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("AllOut/magnolia-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()


################# The End ########################
