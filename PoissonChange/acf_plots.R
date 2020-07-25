set.seed(1)
library("MCMCpack")
library(rep.acf.ccf)
load("mida.rda")

c0 <- 13
d0 <- 1
bars <- 6
m <- 2       #### Number of chains
mc.chain.list <- list()

for (i in 1:m)
  (mc.chain.list[[i]] <- MCMCpoissonChange(mida~ 1, m=bars, c0=c0, d0=d0, marginal.likelihood="none", mcmc = 1e5, burnin = 0, thin = 1, seed = round(1e3*runif(1))))

save(mc.chain.list, file = "Out/two_chains.Rdata")

load(file = "Out/two_chains.Rdata")
lag.max <- 150

###########################################
###########################################
######## ACF and G-ACF ###################
##########################################
###########################################


### ncrop = 1e3

ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)

for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), width = 8, height = 4)
par(mfrow = c(1,2))
plot(local.acf, main = expression(paste("Old ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
plot(global.acf, main = expression(paste("New ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
dev.off()

### ncrop = 1e4


ncrop <- 1e4
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)

for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), width = 8, height = 4)
par(mfrow = c(1,2))
plot(local.acf, main = expression(paste("Old ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
plot(global.acf, main = expression(paste("New ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
dev.off()

### ncrop = 1e5


ncrop <- 1e5
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)

for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), width = 8, height = 4)
par(mfrow = c(1,2))
plot(local.acf, main = expression(paste("Old ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[1]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[1]]$acf), max(true.acf[1,1,]))))
plot(global.acf, main = expression(paste("New ACF for chain 1")), xlim = range(0, lag.max), ylim =  c(0,1))#c(min(min(acf.list[[1]][[2]]$acf), min(true.acf[1,1,])), max(max(acf.list[[1]][[2]]$acf), max(true.acf[1,1,]))))
dev.off()
