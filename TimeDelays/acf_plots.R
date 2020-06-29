library(timedelay)
library(rep.acf.ccf)
source("functions.R")
data(simple)


lcA <- simple[, 1 : 3]
lcB <- simple[, c(1, 4, 5)]
m <- 5
delta.start <- c(10, 30, 50, 70, 90)
micro <- 0
mc.chain.list <- list()

for (i in 1:m){
  mc.chain.list[[i]] <- markov.chain(1e5, delta.start[i], lcA, lcB, micro)$samples
}

########################################
ncrop <- 5e2
########################################

x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

acf <- combined_acf(x, chain = 1, component = 1, lag.max  = 100, type = "correlation")
pdf(file = paste("Out/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
pdf(file = paste("Out/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
pdf(file = paste("Out/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
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
pdf(file = paste("Out/acf_n", ncrop, "_component - 1", ".pdf", sep = ""), height = 3)
par(mfrow = c(1,2))
plot(acf[[1]][[1]], main = expression("Locally centered ACF"))
plot(acf[[1]][[2]], main = expression("Globally centered ACF"))
dev.off()

####################################################

pdf(file = "Out/density _5e2_1e4.pdf")
par(mfrow = c(5,2))
plot(density(mc.chain.list[[1]][1:5e2,1]))
plot(density(mc.chain.list[[1]][1:99999,1]))

plot(density(mc.chain.list[[2]][1:5e2,1]))
plot(density(mc.chain.list[[2]][1:99999,1]))

plot(density(mc.chain.list[[3]][1:5e2,1]))
plot(density(mc.chain.list[[3]][1:99999,1]))

plot(density(mc.chain.list[[4]][1:5e2,1]))
plot(density(mc.chain.list[[4]][1:99999,1]))

plot(density(mc.chain.list[[5]][1:5e2,1]))
plot(density(mc.chain.list[[5]][1:99999,1]))

dev.off()

pdf(file = "Out/trace _5e2_1e4.pdf")
par(mfrow = c(5,2))
plot.ts(mc.chain.list[[1]][1:5e2,1])
plot.ts(mc.chain.list[[1]][1:9999,1])

plot.ts(mc.chain.list[[2]][1:5e2,1])
plot.ts(mc.chain.list[[2]][1:9999,1])

plot.ts(mc.chain.list[[3]][1:5e2,1])
plot.ts(mc.chain.list[[3]][1:9999,1])

plot.ts(mc.chain.list[[4]][1:5e2,1])
plot.ts(mc.chain.list[[4]][1:9999,1])

plot.ts(mc.chain.list[[5]][1:5e2,1])
plot.ts(mc.chain.list[[5]][1:9999,1])

dev.off()

a <- markov.chain(1e6, 30, lcA, lcB, 0)
a <- a$samples
pdf(file = "Out/true_density_plot_n=1e6.pdf")
par(mfrow = c(1,1))
plot(density(a[,1]))
dev.off
