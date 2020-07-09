set.seed(100)
library(timedelay)
library(rep.acf.ccf)
source("functions.R")
data <- read.csv("q0957usno.csv", header = TRUE)


lcA <- data[, 1 : 3]
lcB <- data[, c(1, 4, 5)]
m <- 4
delta.start <- c(300, 500, 700, 900)
micro <- 0
mc.chain.list <- list()

for (i in 1:m){
  a <- markov.chain(lcA, lcB, 1e5+1, delta.start[i], delta.jump = 1, micro = 0, ram = TRUE, burn = 0)
  mc.chain.list[[i]] <- a$samples
  print(a$delta.accept.rate)
}


########################################
ncrop <- 5e2
########################################

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
ncrop <- 1e3
########################################

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
ncrop <- 1e4
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop-1,]
}
acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
pdf(file = paste("Out/acf_ram_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

####################################################
par(mfrow = c(4,2))
for (i in 1:m){
  pdf(file = paste("Out/density_ram_", i, ".pdf", sep = ""), height = 3)
  
  plot(density(mc.chain.list[[i]][1:5e2,1]), xlab = "Time Delay", ylab = "Density", main = "nsim = 500")
  plot(density(mc.chain.list[[i]][1:1e5,1]), xlab = "Time Delay", ylab = "Density", main = "nsim = 100000")
  dev.off()
}


par(mfrow = c(4,2))
for (i in 1:m){
  pdf(file = paste("Out/trace__ram_", i, ".pdf", sep = ""), height = 3)
  
  plot.ts(mc.chain.list[[i]][1:5e2,1], ylab = "Time Delay", main = "nsim = 500")
  plot.ts(mc.chain.list[[i]][1:1e5,1], ylab = "Time Delay", main = "nsim = 100000")
  dev.off()
}

################################################
#### Testing kit
################################################

a <- markov.chain(1e5, 400, lcA, lcB, 0, delta.jump = 100, ram = FALSE)
b <- markov.chain(1e5, 400, lcA, lcB, 0, delta.jump = 100, ram = TRUE)
a <- a$samples
b <- b$samples
plot(density(a[1:1e3,1]))
plot(density(b[1:1e3,1]))
plot.ts(a[,1])
plot.ts(b[,1])

