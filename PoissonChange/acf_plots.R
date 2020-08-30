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

###########################################
###########################################
######## ACF and G-ACF ###################
##########################################
###########################################

component <- 2
### ncrop = 1e3

ncrop <- 1e3
x <- list()

for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("Out/poisson-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

### ncrop = 1e4


ncrop <- 1e4
x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("Out/poisson-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

#########################################
#########################################
########### Trace plots #################
#########################################
#########################################

n <- 5e4
pdf(file = paste("Out/poisson-trace_n", n, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot.ts(mc.chain.list[[1]][1:n,2], col = "steelblue1", xlab = "Time", ylab = "Component-2", main = "")
lines(1:n, mc.chain.list[[2]][1:n, 2], col = "dodgerblue4", xlab = "Time", ylab = "Component-2", main = "")

plot.ts(mc.chain.list[[1]][1:n,3], col = "steelblue1", xlab = "Time", ylab = "Component-3", main = "")
lines(1:n, mc.chain.list[[2]][1:n,3], col = "dodgerblue4", xlab = "Time", ylab = "Component-3", main = "")

dev.off()
