set.seed(1)
source("functions.R")
library(rep.acf.ccf)

#######################################
##### Setting-1
#######################################

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
  chain <- markov.chain(A, B, C, 1e5, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m
component <- 1


####################################
nsim <- 1000
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf", A, B, C, "n", sep = "_"),  nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("A-ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("G-ACF"))

dev.off()



####################################
nsim <- 1e5
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", component = 1 , graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf_n", A, B, C, "n", sep = "_"), nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("A-ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("G-ACF"))
dev.off()

#####################################
#### Setting-2
#####################################

m <- 5
A <- 1
B <- 10
C <- 7
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}

mc.chain.list <- list()
for (i in 1:m){
  chain <- markov.chain(A, B, C, 1e5, start[i,])
  mc.chain.list[[i]] = chain
}
component <- 1


####################################
nsim <- 1000
########################################

x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf", A, B, C, "n", sep = "_"),  nsim, ".pdf", sep = ""), width = 10, height = 5)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("A-ACF"))
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = expression("G-ACF"))

dev.off()
