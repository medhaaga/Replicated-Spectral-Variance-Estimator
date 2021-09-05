set.seed(10)
source("functions.R")
#sourceCpp("MCMC.cpp")
library(multichainACF)


######################################
##### Model parameters ###############
######################################

m <- 5
p <- 100
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(0.99 + seq(0, (p-1))*((0.999 - 0.99)/(p-1)))
dummy <- matrix(1:p^2, nrow = p, ncol = p)
dummy <- qr.Q(qr(dummy))
phi <- dummy %*% phi %*% t(dummy)

target <- target.sigma(phi, omega)
truth <- true.sigma(phi, var = target)
lag.max <- 40

true.acf <- array(0, dim = c(p, p, 2*lag.max + 1))
true.acf[,,lag.max+1] <- target
for (i in 1:lag.max){
  true.acf[,,lag.max + 1 + i] <- phi %*% true.acf[,,lag.max + 1 + i - 1]
  true.acf[,,lag.max + 1 - i] <- true.acf[,,lag.max + 1 - i + 1] %*% t(phi)
}

################ Only for creating Markov chains. Don't run. ##################
start <- matrix(0, nrow = m, ncol = p)  #only depends on C

for(i in 1:max(c(floor(m/2), 1))){
  start[i,] <- p*i*sqrt(diag(target))
  start[m-i+1,] <- -p*i*sqrt(diag(target))
}


mc.chain.list <- list()
start.time = Sys.time()
global.mean <- rep(0,p)
for (i in 1:m){
  #chain <- markov_chain(phi=phi, omega=omega, nsim=5e4, start=start[i,])
  chain <- markov.chain(phi=phi, omega=omega, nsim=5e4, start=start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
  print(colMeans(chain))
}
end.time <- Sys.time()
print(end.time -start.time)


plot(mc.chain.list[[1]][1:1e3,1], mc.chain.list[[1]][1:1e3,p], xlim = c(-1450, 1450), ylim = c(-1450,1450))
points(mc.chain.list[[2]][1:1e3,1], mc.chain.list[[2]][1:1e3,p], col = "red")
points(mc.chain.list[[3]][1:1e3,1], mc.chain.list[[3]][1:1e3,p], col = "green")
points(mc.chain.list[[4]][1:1e3,1], mc.chain.list[[4]][1:1e3,p], col = "blue")
points(mc.chain.list[[5]][1:1e3,1], mc.chain.list[[5]][1:1e3,p], col = "pink")


global.mean <- global.mean/m
save(mc.chain.list, true.acf, file = "Out/var-five_chains.Rdata")
#################################################################################

load(file = "Out/var-five_chains.Rdata")
#######################################################
############### ACF and G-ACF #########################
#######################################################

p <- 100
m <- 5
cp1 <- 1
cp2 <- p
lag.max <- 40

###############################
### ncrop = 1e3
##############################
 
ncrop <- 1e3
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf1 <- globalACF(x, type = "correlation", component = cp1, lag.max = lag.max, chains = c(0), plot = FALSE, avg = FALSE)[[1]]
local.acf1 <- globalACF(x, type = "correlation", component = cp1, mean = "local", lag.max = lag.max, chains = c(0), plot = FALSE, avg = FALSE)[[1]]
global.acf2 <- globalACF(x, type = "correlation", component = cp2, lag.max = lag.max, chains = c(0), plot = FALSE, avg = FALSE)[[1]]
local.acf2 <- globalACF(x, type = "correlation", component = cp2, mean = "local", lag.max = lag.max, chains = c(0), plot = FALSE, avg = FALSE)[[1]]

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf1, main = expression("Locally centered ACF"), col = "blue", type = 'l', lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp1,cp1,]/true.acf[cp1,cp1,lag.max + 1], col = adjustcolor("blue", alpha=.4), lwd=2, lty=2)
lines(seq(0, lag.max), local.acf2[[1]], col = "green3", lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp2,cp2,]/true.acf[cp2,cp2,lag.max + 1], col = adjustcolor("green3", alpha=.4), lwd=2, lty=2)
legend("bottomleft", legend = c("ACF, Comp=1", "True ACF, Comp=1", "ACF, Comp=100", "True ACF, Comp=100"), cex=.7, col = c("blue", "blue", "green3", "green3"), lty=rep(c(1,2),2), lwd=2)

plot(global.acf1, main = expression("Globally centered ACF"), col = "blue", type = 'l', lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp1,cp1,]/true.acf[cp1,cp1,lag.max + 1], col =  adjustcolor("blue", alpha=.5), lwd=2, lty=2)
lines(seq(0, lag.max), global.acf2[[1]], col = "green3", lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp2,cp2,]/true.acf[cp2,cp2,lag.max + 1], col = adjustcolor("green3", alpha=.4), lwd=2, lty=2)
legend("bottomleft", legend = c("G-ACF, Comp=1", "True ACF, Comp=1", "G-ACF, Comp=100", "True ACF, Comp=100"), cex=.7, col = c("blue", "blue", "green3", "green3"), lty=rep(c(1,2),2), lwd=2)

dev.off()

###############################
### ncrop = 5e4
##############################

ncrop <- 5e4
x <- list()
for (i in 1:m){
  x[[i]] <- mc.chain.list[[i]][1:ncrop,]
}

global.acf1 <- globalACF(x, type = "correlation", component = cp1, lag.max = lag.max, chains = c(1), plot = FALSE, avg = FALSE)[[1]]
local.acf1 <- globalACF(x, type = "correlation", component = cp1, mean = "local", lag.max = lag.max, chains = c(1), plot = FALSE, avg = FALSE)[[1]]
global.acf2 <- globalACF(x, type = "correlation", component = cp2, lag.max = lag.max, chains = c(1), plot = FALSE, avg = FALSE)[[1]]
local.acf2 <- globalACF(x, type = "correlation", component = cp2, mean = "local", lag.max = lag.max, chains = c(1), plot = FALSE, avg = FALSE)[[1]]

pdf(file = paste("Out/var-acf_n", ncrop, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,2))
plot(local.acf1, main = expression("Locally centered ACF"), col = "blue", type = 'l', lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp1,cp1,]/true.acf[cp1,cp1,lag.max + 1], col = adjustcolor("blue", alpha=.4), lwd=2, lty=2)
lines(seq(0, lag.max), local.acf2[[1]], col = "green3", lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp2,cp2,]/true.acf[cp2,cp2,lag.max + 1], col = adjustcolor("green3", alpha=.4), lwd=2, lty=2)
legend("bottomleft", legend = c("ACF, Comp=1", "True ACF, Comp=1", "ACF, Comp=100", "True ACF, Comp=100"), cex=.7, col = c("blue", "blue", "green3", "green3"), lty=rep(c(1,2),2), lwd=2)

plot(global.acf1, main = expression("Globally centered ACF"), col = "blue", type = 'l', lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp1,cp1,]/true.acf[cp1,cp1,lag.max + 1], col =  adjustcolor("blue", alpha=.5), lwd=2, lty=2)
lines(seq(0, lag.max), global.acf2[[1]], col = "green3", lwd=2)
lines(seq(-lag.max, lag.max), true.acf[cp2,cp2,]/true.acf[cp2,cp2,lag.max + 1], col = adjustcolor("green3", alpha=.4), lwd=2, lty=2)
legend("bottomleft", legend = c("G-ACF, Comp=1", "True ACF, Comp=1", "G-ACF, Comp=100", "True ACF, Comp=100"), cex=.7, col = c("blue", "blue", "green3", "green3"), lty=rep(c(1,2),2), lwd=2)

dev.off()



