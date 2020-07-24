set.seed(1)
source("functions.R")
library(rep.acf.ccf)


m <- 5
A <- 1
B <- 3
C <- 8
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, 5e4, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m
lag.max <- 100
component <- 1


####################################
nsim <- 1000
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/acf_n", sep = "_"), nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Locally-centered ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Globally-centered ACF"))

dev.off()



####################################
nsim <- 5000
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/acf_n", sep = "_"), nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Locally-centered ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Globally-centered ACF"))

dev.off()



####################################
nsim <- 1e4
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1 , graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/acf_n", sep = "_"), nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Locally-centered ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Globally-centered ACF"))

dev.off()

####################################
nsim <- 5e4
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/acf_n", sep = "_"), nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Locally-centered ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("Globally-centered ACF"))

dev.off()


####################################
nsim <- 1e3
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]

pdf(file = paste("Out/_1_3_8_/globalACF_all_chains_n", nsim, ".pdf", sep = ""))
globalACF(x, type = "correlation", lag.max = lag.max, chains = 0, component = 1)
dev.off()

####################################
nsim <- 1e4
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]

pdf(file = paste("Out/_1_3_8_/globalACF_all_chains_n", nsim, ".pdf", sep = ""))
globalACF(x, type = "correlation", lag.max = lag.max, chains = 0, component = 1)
dev.off()

####################################
nsim <- 5e4
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]

pdf(file = paste("Out/_1_3_8_/globalACF_all_chains_n", nsim, ".pdf", sep = ""))
globalACF(x, type = "correlation", lag.max = lag.max, chains = 0, component = 1)
dev.off()
