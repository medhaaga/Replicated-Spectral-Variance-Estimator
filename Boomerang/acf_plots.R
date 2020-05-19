set.seed(100)
source("functions.R")
library(rep.acf.ccf)


m <- 5
nsim <- 500
A <- 2
B <- 9
C <- 7
start <- matrix(0, nrow = m, ncol = 2)  #only depends on C


for(i in 1:floor(m/2)){
  start[i,] <- c(0, C*(2^(2-i)))
  start[m-i+1,] <- c(C*(2^(2-i)), 0)
}
mc.chain.list <- list()
global.mean <- c(0,0)
for (i in 1:m){
  chain <- markov.chain(A, B, C, nsim, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
}
global.mean <- global.mean/m

acf.list <- combined_acf(mc.chain.list, center = global.mean, chain = 1, component = c(1,2), lag.max = 150, type = "correlation")
ccf.list <- combined_ccf(mc.chain.list, center = global.mean, chain = 1, component = c(1,2), lag.max = 150, type = "correlation")

pdf(file = paste(paste("Out/", A, B, C, "/acf, n=", sep = "_"), nsim, ".pdf", sep = ""), title = paste("Locally vs. globally centered ACF plot, nsim = ", nsim))
par(mfrow = c(3,2))
for (i in 1:2){
  plot(acf.list[[i]][[1]], main = paste("Locally centered ACF plot for component -", i ))
  plot(acf.list[[i]][[2]], main = paste("Globally centered ACF plot for component -", i ))
}
plot(ccf.list[[1]], main = "Locally centered CCF plot")
plot(ccf.list[[2]], main = "Globally centered CCF plot")
dev.off()
