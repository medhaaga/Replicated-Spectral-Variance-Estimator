set.seed(10)
source("functions.R")
library(multichainACF)


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

phi <- diag(c(.999, .001))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
lag.max <- 40

true.acf <- array(0, dim = c(p, p, 2*lag.max + 1))
true.acf[,,lag.max+1] <- target
for (i in 1:lag.max){
  true.acf[,,lag.max + 1 + i] <- phi %*% true.acf[,,lag.max + 1 + i - 1]
  true.acf[,,lag.max + 1 - i] <- true.acf[,,lag.max + 1 - i + 1] %*% t(phi)
}

################ Only for creating Markov chains. Don't run. ##################
start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in 1:floor(m/2)){
  start[i,] <- i*sqrt(diag(target))
  start[m-i+1,] <- -i*sqrt(diag(target))
}

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  chain <- markov.chain(phi, omega, 5e4, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
  print(colMeans(chain))
}
global.mean <- global.mean/m
save(mc.chain.list, true.acf, file = "Out/var-five_chains.Rdata")
#################################################################################

load(file = "Out/var-five_chains.Rdata")
#######################################################
############### ACF and G-ACF #########################
#######################################################

component <- 1
lag.max <- 40

###############################
### ncrop = 1e3
##############################

ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", component = 1, lag.max = lag.max, chains = c(2), plot = FALSE, avg = FALSE)[[1]]
local.acf <- globalACF(x, type = "correlation", component = 1, mean = "local", lag.max = lag.max, chains = c(2), plot = FALSE, avg = FALSE)[[1]]

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
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

global.acf <- globalACF(x, type = "correlation", component = 1, lag.max = lag.max, chains = c(2), plot = FALSE, avg = FALSE)[[1]]
local.acf <- globalACF(x, type = "correlation", component = 1, mean = "local", lag.max = lag.max, chains = c(2), plot = FALSE, avg = FALSE)[[1]]

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
plot(global.acf, main = expression("Globally centered ACF"))
lines(seq(-lag.max, lag.max), true.acf[1,1,]/true.acf[1,1,lag.max + 1], col = "red")
dev.off()
