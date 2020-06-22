set.seed(130290)
library(mcmcse)
library(mcmc)

source("spatio_temp_function.R")
load("spatio_temp_truth_2e6by1000")
truth <- colMeans(spatio_truth)
data(logit)
level <- .90

## Fitting model with intercept
N <- 500
p <- 2

starting1 <- list("beta"=rep(30,N.t*p), "phi"=rep(3/(0.07*max.d), N.t),
                  "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                  "sigma.eta"=diag(rep(0.01, p)))

starting2 <- list("beta"=rep(-30,N.t*p), "phi"=rep(3/(0.89*max.d), N.t),
                  "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                  "sigma.eta"=diag(rep(0.01, p)))

chain1 <- spatio_temp(N = N, verbose = FALSE, starting1)
chain2 <- spatio_temp(N = N, verbose = FALSE, starting2)

chains <- list(chain1, chain2)

plot(chain1[,1], chain1[,11], col = "blue", xlim = range(chain1[,1], 
    chain1[,11], chain2[,1], chain2[,11]), ylim = range(chain1[,1], 
                                                         chain1[,11], chain2[,1], chain2[,11]))
points(chain2[,1], chain2[,11], col = "green", xlim = range(chain1[,1], 
                                                           chain1[,11], chain2[,1], chain2[,11]), ylim = range(chain1[,1], 
                                                                                                             chain1[,11], chain2[,1], chain2[,11]))

#

##############################
ncrop <- 5e2
##############################

acf <- combined_acf(list(chains[[1]][1:ncrop,], chains[[2]][1:ncrop,]), chain = 1, component = (1:26), 
                    lag.max = 100, type = "correlation")

for (j in 1:26)
{
  pdf(file = paste("Out/acf_n", ncrop, "_component", j, ".pdf", sep = ""), height = 3, width = 4)
  plot(acf[[j]][[1]], main = expression("Autocorrelation function"), type = "l")
  lines(seq(0:lag.max), acf[[j]][[2]]$acf, main = expression("Globally centered ACF"), col = "red")
  legend("topright", legend=c("ACF", "R-ACF"),col=c("black", "red"), lty=1, cex=.5)
  dev.off()
}


##############################
ncrop <- 1e4
##############################

acf <- combined_acf(list(chains[[1]][1:ncrop,], chains[[2]][1:ncrop,]), chain = 1, component = (1:26), 
                    lag.max = 100, type = "correlation")

for (j in 1:26)
{
  pdf(file = paste("Out/acf_n", ncrop, "_component", j, ".pdf", sep = ""), height = 3, width = 4)
  plot(acf[[j]][[1]], main = expression("Autocorrelation function"), type = "l")
  lines(seq(0:lag.max), acf[[j]][[2]]$acf, main = expression("Globally centered ACF"), col = "red")
  legend("topright", legend=c("ACF", "R-ACF"),col=c("black", "red"), lty=1, cex=.5)
  dev.off()
}