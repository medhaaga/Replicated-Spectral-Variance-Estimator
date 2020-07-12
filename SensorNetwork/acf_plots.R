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

nsim = 1e5

start1 <- runif(8, min= -0.3, max = 0)
aux1 <- runif(8, min = -0.3, max= 0)
system.time(res.ram1 <- MHwG.RAM(start1, aux1, jump.scale = j.scale,
                                 Ob, Os, Xb, Xs, Yb, Ys,
                                 n.sample = nsim, n.burn = 0))

# Chain-2 with 5e3 samples
start2 <- runif(8, min = 0.7, max = 1.1)
aux2 <- runif(8, min = 0.7, max= 1.1)

system.time(res.ram2 <- MHwG.RAM(start2, aux2, jump.scale = j.scale,
                                    Ob, Os, Xb, Xs, Yb, Ys,
                                    n.sample = nsim, n.burn = 0))

save(res.ram1, res.ram2, file = "Out/two_chains.Rdata")
colMeans(mc.chain.list[[1]])
colMeans(mc.chain.list[[2]])

mc.chain.list <- list(res.ram1$x, res.ram2$x)

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

########################################
ncrop <- 1e3
########################################
m <- 2
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 5e3
########################################
m <- 2
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component =7, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 7", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 1e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component =7, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 7", ".pdf", sep = ""), height = 3)

par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 5e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop-1,]
}
acf <- combined_acf(x, chain = 1, component =7, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 7", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 1e5-1
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component =7, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 7", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()



# Figure 6

par(mfrow = c(4, 3), font = 2, font.lab = 2, font.axis = 2, cex = 1.2,
    mai = c(0.7, 0.9, 0.7, 0.3), mgp = c(2.5, 0.5, 0), las = 1)


hist(res.ram$x[, 1], 50, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 11), xlim = c(-0.2, 0.75))
title(expression(bold(paste("RAM: ", x[11]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
index <- which(!is.na(res))
lines(density(res.ram.den$x[, 1]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")),
                expression(bold("True location"))), lwd = 3,
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.5748, lty = 2, lwd = 2)


hist(res.ram$x[, 3], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 3), xlim = c(-0.15, 1.3))
title(expression(bold(paste("RAM: ", x[21]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.8, cex = 1.2, las = 0)
legend("top", c(expression(bold("Marginal dist.")),
                expression(bold("True location"))), lwd = 3,
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
lines(density(res.ram.den$x[, 3]), lwd = 1.5)
abline(v = 0.1, lty = 2, lwd = 2)


hist(res.ram$x[, 5], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 15), xlim = c(0.15, 0.85))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
index <- which(!is.na(res))
title(expression(bold(paste("RAM: ", x[31]))))
lines(density(res.ram.den$x[, 5]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")),
                expression(bold("True location"))), lwd = 3,
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.2578, lty = 2, lwd = 2)


hist(res.ram$x[, 7], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 2.3), xlim = c(-1.1, 1.8))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[41]))))
lines(density(res.ram.den$x[, 7]), lwd = 1.5)
legend("topleft", c(expression(bold("Marginal dist.")),
                    expression(bold("True location"))), lwd = 3,
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.85, lty = 2, lwd = 2)

