set.seed(1)
source("functions.R")
library(multichainACF)

#######################################
##### Setting-1
#######################################

m <- 5
A <- 1
B <- 3
C <- 8

############# Run only for simulating Markov chains. ############
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
save(mc.chain.list, file = "Out/_1_3_8_/boom-five_chains_1_3_8.Rdata")
#####################################################################

load(file = "Out/_1_3_8_/boom-five_chains_1_3_8.Rdata")
component <- 1
lag.max <- 50

####################################
nsim1 <- 1000
nsim2 <- 1e5
########################################

x <- list()
y <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:nsim1,]
  y[[i]] <- mc.chain.list[[i]][1:nsim2,]
}

global.acf1 <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, plot = FALSE)$"avgACF"
global.acf2 <- globalACF(y, type = "correlation", lag.max = lag.max, component = 1, plot = FALSE)$"avgACF"
local.acf1 <- globalACF(x, type = "correlation", lag.max = lag.max, mean = "local", component = 1, plot = FALSE)$"avgACF"
local.acf2 <- globalACF(y, type = "correlation", lag.max = lag.max, mean = "local", component = 1, plot = FALSE)$"avgACF"

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf", A, B, C, sep = "_"), ".pdf", sep = ""), width = 10, height = 4)
par(mfrow = c(1,2))
plot(local.acf1, xlab = "Lag", ylab = "Autocorrelation", main = "")
lines(local.acf2$acf, col = "steelblue1", lwd=2)
plot(global.acf1, xlab = "Lag", ylab = "Autocorrelation", main = "")
lines(global.acf2$acf, col = "steelblue1", lwd=2)
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

save(mc.chain.list, file = "Out/_1_10_7_/boom-five_chains_1_10_7.Rdata")

load(file = "Out/_1_10_7_/boom-five_chains_1_10_7.Rdata")
####################################
nsim <- 1000
########################################


x <- list()
for (i in 1:m)
  x[[i]] <- mc.chain.list[[i]][1:nsim,]


global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, plot = FALSE)$"avgACF"
local.acf <- globalACF(x, type = "correlation", lag.max = lag.max, mean = "local", component = 1, plot = FALSE)$"avgACF"

pdf(file = paste(paste("Out/", A, B, C, "/boom-acf", A, B, C, "n", sep = "_"),  nsim, ".pdf", sep = ""), width = 10, height = 4)
par(mfrow = c(1,2))
plot(local.acf, xlab = "Lag", ylab = "Autocorrelation", main = "")
plot(global.acf, xlab = "Lag", ylab = "Autocorrelation", main = "")
dev.off()
