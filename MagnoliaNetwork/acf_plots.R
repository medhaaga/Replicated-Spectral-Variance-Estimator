
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


########################################
ncrop <- 1e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

acf <- combined_acf(x, chain = 1, component = 5, lag.max  = 50, type = "correlation")

pdf(file = paste("Out/acf_n", ncrop, "_component - 5", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("ACF for race"))
plot(acf[[1]][[2]], main = expression("R-ACF for race"))
dev.off()


########################################
ncrop <- 1e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component = 5, lag.max  = 50, type = "correlation")

pdf(file = paste("Out/acf_n", ncrop, "_component - 5", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 1e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop-1,]
}
acf <- combined_acf(x, chain = 1, component = 5, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_n", ncrop, "_component - 5", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("ACF for race"))
plot(acf[[1]][[2]], main = expression("R-ACF for race"))
dev.off()


