set.seed(10)
source("functions.R")
library(multichainACF)
library("mvtnorm")

######################################
##### Model parameters ###############
######################################

m <- 101
p <- 2
omega <- diag(p)
for (i in 1:(p-1)){
  for (j in 1:(p-i)){
    omega[j, j+i] <- .9^i
    omega[j+i, j] <- .9^i
  }
}

phi <- diag(c(.999, .001))
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
for(i in 1:m)
  start[i,] <- rmvnorm(n=1, mean = rep(0,p), sigma = target)
# for(i in 1:floor(m/2)){
#   start[(2*i),] <- i*.5*sqrt(diag(target))
#   start[(2*i +1),] <- -i*.5*sqrt(diag(target))
# }

mc.chain.list <- list()
global.mean <- rep(0,p)
for (i in 1:m){
  print(i)
  print(paste("Starting at: ", start[i,]))
  chain <- markov.chain(phi, omega, 1e4, start[i,])
  global.mean <- global.mean + colMeans(chain)
  mc.chain.list[[i]] = chain
  print(paste("Local mean: ", colMeans(chain)))
}
global.mean <- global.mean/m


save(mc.chain.list, true.acf, file = "Out/var-100_chains.Rdata")
#################################################################################

load(file = "Out/var-100_chains.Rdata")
#######################################################
############### ACF and G-ACF #########################
#######################################################

component <- 1

lag1 <- 1
lag2 <- lag.max

lacf_twolags <- matrix(0, nrow = 2, ncol = floor(m/2))
gacf_twolags <- matrix(0, nrow = 2, ncol = floor(m/2))

nsims <- c(1e3, 1e4)

for (n in 1:2)
{
  ncrop <- nsims[n]
  
  for (t in (1:floor(m/2)))
  {
    print(t)
    nchains <- 2*t
    
    
    x <- list()
    for (i in 1:nchains){
      x[[i]] <- mc.chain.list[[i]][1:ncrop,]
    }
    
    global.acf <- globalACF(x, type = "correlation", component = 1, lag.max = lag.max, chains = 0, plot = FALSE, avg = TRUE)
    local.acf <- globalACF(x, type = "correlation", component = 1, mean = "local", lag.max = lag.max, chains = 0, plot = FALSE, avg = TRUE)
    
    lacf_twolags[1,t] <- local.acf$avgACF[[1]][1+lag1]
    lacf_twolags[2,t] <- local.acf$avgACF[[1]][1+lag2]
    gacf_twolags[1,t] <- global.acf$avgACF[[1]][1+lag1]
    gacf_twolags[2,t] <- global.acf$avgACF[[1]][1+lag2]
    
  }
  
  pdf(file = paste("Out/ACFvsm_", ncrop, ".pdf", sep = ""), height=5, width=5)
  plot(seq(2, 2*floor(m/2), 2), lacf_twolags[1,], type = "l", ylim = c(0,1), col = "blue", lwd=2, ylab = "ÄCF", xlab = "m", main = paste("n =", ncrop))
  lines(seq(2, 2*floor(m/2), 2), lacf_twolags[2,], col = "blue", lwd=2, lty=2)
  lines(seq(2, 2*floor(m/2), 2), gacf_twolags[1,], col = "orange", lwd=2)
  lines(seq(2, 2*floor(m/2), 2), gacf_twolags[2,], col = "orange", lwd=2, lty=2)
  abline(h = true.acf[1,1,lag1+lag.max+1]/true.acf[1,1,lag.max+1], col = "red", lwd=1)
  abline(h = true.acf[1,1,lag2+lag.max+1]/true.acf[1,1,lag.max+1], col = "red", lwd=1, lty=2)
  legend("bottomright", legend = c("A-ACF, lag=1", "G-ACF, lag=1", "Truth, lag=1", "A-ACF, lag=40", "G-ACF, lag=40", "Truth, lag=40"), col = rep(c("blue", "orange", "red"), 2), lty=c(1,1,1,2,2,2), cex=.7)
  dev.off()
}

