set.seed(1)
library(rep.acf.ccf)
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
mc.chain.list <- list()

for (k in 1:m){
  temp <- MHwG.RAM(start[k,], aux[k,], jump.scale = j.scale,Ob, Os, Xb, Xs, Yb, Ys, n.sample = 1e4, n.burn = 0)
  mc.chain.list[[k]] <- temp$x
  print(colMeans(temp$accept))
}

save(mc.chain.list, file = "Out/five_chains.Rdata")

###############################################
###############################################
# Scatter plots
###############################################
###############################################
nsim1 <- 5e3-1
nsim2 <- 1e5-1
min <- min(min(apply(mc.chain.list[[1]], 2, min)), min(apply(mc.chain.list[[2]], 2, min)))
max <- max(max(apply(mc.chain.list[[1]], 2, max)), max(apply(mc.chain.list[[2]], 2, max)))
pdf(file = "Out/scatter_plots.pdf")
par(mfrow = c(4, 2))

# X_1 for nsim1
plot(mc.chain.list[[1]][1:nsim1, c(1, 2)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,1], mc.chain.list[[2]][1:nsim1,1]),
     ylim = range(mc.chain.list[[1]][1:nsim1,2], mc.chain.list[[2]][1:nsim1,2]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(1,2)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,1], mc.chain.list[[2]][1:nsim1,1]),
       ylim = range(mc.chain.list[[1]][1:nsim1,2], mc.chain.list[[2]][1:nsim1,2]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)

#X_1 for nsim2
plot(mc.chain.list[[1]][1:nsim2, c(1, 2)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim2,1], mc.chain.list[[2]][1:nsim2,1]),
     ylim = range(mc.chain.list[[1]][1:nsim2,2], mc.chain.list[[2]][1:nsim2,2]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim2,c(1,2)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim2,1], mc.chain.list[[2]][1:nsim2,1]),
       ylim = range(mc.chain.list[[1]][1:nsim2,2], mc.chain.list[[2]][1:nsim2,2]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)

#X_2 for nsim1
plot(mc.chain.list[[1]][1:nsim1, c(3, 4)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,3], mc.chain.list[[2]][1:nsim1,3]),
     ylim = range(mc.chain.list[[1]][1:nsim1,4], mc.chain.list[[2]][1:nsim1,4]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(3,4)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,3], mc.chain.list[[2]][1:nsim1,3]),
       ylim = range(mc.chain.list[[1]][1:nsim1,4], mc.chain.list[[2]][1:nsim1,4]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)

#X_2 for nsim2

plot(mc.chain.list[[1]][1:nsim1, c(3, 4)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,3], mc.chain.list[[2]][1:nsim1,3]),
     ylim = range(mc.chain.list[[1]][1:nsim1,4], mc.chain.list[[2]][1:nsim1,4]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(3,4)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,3], mc.chain.list[[2]][1:nsim1,3]),
       ylim = range(mc.chain.list[[1]][1:nsim1,4], mc.chain.list[[2]][1:nsim1,4]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)

#X_3 for nsim1

plot(mc.chain.list[[1]][1:nsim1, c(5, 6)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,5], mc.chain.list[[2]][1:nsim1,5]),
     ylim = range(mc.chain.list[[1]][1:nsim1,6], mc.chain.list[[2]][1:nsim1,6]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(5,6)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,5], mc.chain.list[[2]][1:nsim1,5]),
       ylim = range(mc.chain.list[[1]][1:nsim1,6], mc.chain.list[[2]][1:nsim1,6]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)


#X_3 for nsim2

plot(mc.chain.list[[1]][1:nsim2, c(5, 6)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim2,5], mc.chain.list[[2]][1:nsim2,5]),
     ylim = range(mc.chain.list[[1]][1:nsim2,6], mc.chain.list[[2]][1:nsim2,6]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim2,c(5,6)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim2,5], mc.chain.list[[2]][1:nsim2,5]),
       ylim = range(mc.chain.list[[1]][1:nsim2,6], mc.chain.list[[2]][1:nsim2,6]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)

#X_4 for nsim1

plot(mc.chain.list[[1]][1:nsim1, c(7, 8)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,7], mc.chain.list[[2]][1:nsim1,7]),
     ylim = range(mc.chain.list[[1]][1:nsim1,8], mc.chain.list[[2]][1:nsim1,8]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(7,8)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,7], mc.chain.list[[2]][1:nsim1,7]),
       ylim = range(mc.chain.list[[1]][1:nsim1,8], mc.chain.list[[2]][1:nsim1,8]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)

#X_4 for nsim2

plot(mc.chain.list[[1]][1:nsim1, c(7, 8)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,7], mc.chain.list[[2]][1:nsim1,7]),
     ylim = range(mc.chain.list[[1]][1:nsim1,8], mc.chain.list[[2]][1:nsim1,8]), xlab = "", ylab = "", main = "")
points(mc.chain.list[[2]][1:nsim1,c(7,8)], pch = 1, xlim = range(mc.chain.list[[1]][1:nsim1,7], mc.chain.list[[2]][1:nsim1,7]),
       ylim = range(mc.chain.list[[1]][1:nsim1,8], mc.chain.list[[2]][1:nsim1,8]), xlab = "", ylab = "", main = "", col = "red")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)

dev.off()
#####################################################
#####################################################
### ACF plots
#####################################################
#####################################################
component <- 1
########################################
ncrop <- 1e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 5e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 1e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 5e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 1e5
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf <- globalACF(x, type = "correlation", lag.max = lag.max, component = 1, graph = FALSE)$'G-ACF'
local.acf <- acf(x[[1]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)
for (i in 2:m)
  local.acf$acf <- acf(x[[i]][, component], type = "correlation", lag.max = lag.max, plot = FALSE)$acf
local.acf$acf <- local.acf$acf/m

pdf(file = paste("Out/acf_n", ncrop, ".pdf", sep = ""), height = 5, width = 10)
par(mfrow = c(1,2))
plot(local.acf, main = expression("Locally centered ACF"))
plot(global.acf, main = expression("Globally centered ACF"))
dev.off()
