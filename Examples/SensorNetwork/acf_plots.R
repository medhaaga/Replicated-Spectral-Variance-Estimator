set.seed(1)
library(multichainACF)
source("functions.R")
######## Data

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors.
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)

# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546,
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0,
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)


j.scale <- rep(1.08, 4)    # jump scale

# Chain-1 with 5e3 samples
m <- 5
start1 <- c(-0.1, 0.5, -0.1, -0.2, 0.1, 0.1, -0.5, -0.5)
aux1 <- runif(n=8, min=min(start1), max=max(start1))
start2 <- c(0.0, 0.6, 0.1, 0.1, 0.2, 0.2, 1.0, 0.0)
aux2 <- runif(n=8, min=min(start2), max=max(start2))
start3 <- c(0.2, 0.7, 0.5, 0.4, 0.5, 0.3, 0.5, 0.5)
aux3 <- runif(n=8, min=min(start3), max=max(start3))
start4 <- c(0.4, 0.8, 0.8, 0.6, 0.7, 0.4, 1.0, 1.0)
aux4 <- runif(n=8, min=min(start4), max=max(start4))
start5 <- c(0.7, 1.0, 1.2, 0.9, 0.9, 0.5, 1.5, 1.5)
aux5 <- runif(n=8, min=min(start5), max=max(start5))

start <- rbind(start1, start2, start3, start4, start5)
aux <- rbind(aux1, aux2, aux3, aux4, aux5)
j.scale <- rep(0.5, 4)

############# Run only for simulating marov chains #####################
mc.chain.list <- list()

for (k in 1:m){
  temp <- MHwG.RAM(start[k,], aux[k,], jump.scale = j.scale,Ob, Os, Xb, Xs, Yb, Ys, n.sample = 1e5+1, n.burn = 0)
  mc.chain.list[[k]] <- temp$x
  print(colMeans(temp$accept))
}

save(mc.chain.list, file = "Out/five_chains.Rdata")
######################################################################

load(file = "Out/five_chains.Rdata")

#####################################################
#####################################################
### ACF plots
#####################################################
#####################################################

component <- 1
lag.max <- 40

########################################
ncrop <- 5e3
########################################
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}


pdf(file = "Out/sensor-acf_n5e3.pdf", height = 4, width = 4)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "local", type = "correlation", col = "darkorange")
par(new = TRUE)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "global", type = "correlation", col = "royalblue")
dev.off()


########################################
ncrop <- 5e4
########################################
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}


pdf(file = "Out/sensor-acf_n5e4.pdf", height = 4, width = 4)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "local", type = "correlation", col = "darkorange")
par(new = TRUE)
globalACF(x, chains = 0, component = component, lag.max = lag.max, mean = "global", type = "correlation", col = "royalblue")
dev.off()


