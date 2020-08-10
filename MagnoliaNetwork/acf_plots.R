
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
  set.seed(1)
  m <- 2
  start <- c(386, 244)
  markovChains = magnolia.mhrw.func(nrep = 1e3, chains = m, start = start)
  mc.chain.list <- list()
  for(i in 1:m){
    mc.chain.list[[i]] <- MHRWestimation(markovChains[[i]])
    mc.chain.list[[i]] = mc.chain.list[[i]][,c(1,7,5,8,9)]
  }

  component = 3

########################################
ncrop <- 1e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("Out/magnolia-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()

#######################################
ncrop <- 1e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- as.matrix(mc.chain.list[[i]][1:ncrop,])
}

pdf(file = paste("Out/magnolia-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
globalACF(x, chains=0, component = component, mean = "local", type = "correlation", leg = FALSE, col = "darkorange")
globalACF(x, chains=0, component = component, mean = "global", type = "correlation", leg = FALSE, col = "royalblue")
dev.off()
