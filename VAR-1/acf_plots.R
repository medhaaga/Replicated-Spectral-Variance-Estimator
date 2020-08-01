set.seed(10)
source("functions.R")
library(rep.acf.ccf)


######################################
##### Model parameters ###############
######################################

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

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)


true.acf <- array(0, dim = c(p, p, 2*lag.max + 1))
true.acf[,,lag.max+1] <- target
for (i in 1:lag.max){
  true.acf[,,lag.max + 1 + i] <- phi %*% true.acf[,,lag.max + 1 + i - 1]
  true.acf[,,lag.max + 1 - i] <- true.acf[,,lag.max + 1 - i + 1] %*% t(phi)
}

start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- 0.5*i*sqrt(diag(target))
  start[m-i+1,] <- -0.5*i*sqrt(diag(target))
}

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, 1e5, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
  print(colMeans(chain))
}
global.mean <- global.mean/m
save(mc.chain.list, file = "Out/five_chains.Rdata")

load(file = "Out/five_chains.Rdata")
#######################################################
############### ACF and G-ACF #########################
#######################################################

component <- 1
###############################
### ncrop = 1e3
##############################

ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()

###############################
### ncrop = 1e4
##############################

ncrop <- 1e4
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)

par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()

###############################
### ncrop = 1e5
##############################

ncrop <- 1e5
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()
