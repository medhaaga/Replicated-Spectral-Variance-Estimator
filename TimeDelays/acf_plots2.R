library(timedelay)
library(rep.acf.ccf)
source("functions.R")
load("dataset")


lcA <- data[, 1 : 3]
lcB <- data[, c(1, 4, 5)]
m <- 5
delta.start <- c(10, 30, 50, 70, 90)
micro <- 0
mc.chain.list <- list()

for (i in 1:m){
  mc.chain.list[[i]] <- markov.chain(1e5+1, delta.start[i], lcA, lcB, micro)$samples
}

########################################
ncrop <- 5e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 1e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

########################################
ncrop <- 5e3
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()


########################################
ncrop <- 1e5
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop-1,]
}
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out2/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

####################################################

for (i in 1:m){
  pdf(file = paste("Out2/density_", i, ".pdf", sep = ""), height = 3)
  par(mfrow = c(1,2))
  plot(density(mc.chain.list[[i]][1:5e2,1]), xlab = "Time Delay", ylab = "Density", main = "nsim = 500")
  plot(density(mc.chain.list[[i]][1:1e5,1]), xlab = "Time Delay", ylab = "Density", main = "nsim = 100000")
  dev.off()
}



for (i in 1:m){
  pdf(file = paste("Out2/trace_", i, ".pdf", sep = ""), height = 3)
  par(mfrow = c(1,2))
  plot.ts(mc.chain.list[[i]][1:5e2,1], ylab = "Time Delay", main = "nsim = 500")
  plot.ts(mc.chain.list[[i]][1:1e5,1], ylab = "Time Delay", main = "nsim = 100000")
  dev.off()
}

a <- markov.chain(1e6, 30, lcA, lcB, 0)
a <- a$samples
pdf(file = "Out2/true_density_plot_n=1e6.pdf")
par(mfrow = c(1,1))
plot(density(a[,1]))
dev.off
