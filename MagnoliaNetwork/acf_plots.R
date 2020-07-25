set.seed(1)
library(rep.acf.ccf)
source("functions.R")
library(ergm)

data(faux.magnolia.high)
magnolia = faux.magnolia.high

#Ensure well connected component only
notwellconnected = which(component.largest(magnolia, connected="weak")==FALSE)
delete.vertices(magnolia, notwellconnected)
dp = dataPrepHighSchool(magnolia)

#get an assortativy measure
gmagnolia = intergraph::asIgraph(magnolia)
igraph::assortativity_degree(gmagnolia)

#m chains length 1e5
  set.seed(39)
  m <- 2
  markovChains = magnolia.mhrw.func(nrep = 1e4, chains = m, start = c(279, 360))
  mc.chain.list <- list()
  for(i in 1:m){
    mc.chain.list[[i]] <- MHRWestimation(markovChains[[i]])
    mc.chain.list[[i]] = mc.chain.list[[i]][,c(1,7,5,8,9)]
  }

lag.max <- 100
component <- 3


########################################
ncrop <- 1e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = component, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 5e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = component, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()

#######################################
ncrop <- 1e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = component, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()

#######################################
ncrop <- 1e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = component, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- local.acf$acf + acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()
