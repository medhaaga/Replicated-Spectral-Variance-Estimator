
library("MCMCpack")
load("mida.rda")

out1 <- MCMCpoissonChange(mida~ 1, m=6, c0=13, d0=1,  
		marginal.likelihood="none", mcmc = 1e4, burnin = 0, thin = 1, seed = round(1e3*runif(1)))

out2 <- MCMCpoissonChange(mida~ 1, m=6, c0=13, d0=1,  
  marginal.likelihood="none", mcmc = 1e4, burnin = 0, seed = round(1e3*runif(1)))

par(mfrow = c(2,7))
for(i in 1:7)
{
  plot.ts(out1[,i])
  plot.ts(out2[,i])
}

par(mfrow = c(2,7))
for(i in 1:7)
{
  acf(out1[,i])
  acf(out2[,i])
}


cbind(colMeans(out1), colMeans(out2))

par(mfrow = c(2,7))
for(i in 1:7)
{
  plot(density(out1[,i]))
plot(density(out2[,i]))
}

ncrop <- 1e3
foo <- list( matrix(out1[1:ncrop, ], nrow = ncrop, ncol = 7) , matrix(out2[1:ncrop, ], nrow = ncrop, ncol = 7))
acfss <- combined_acf(foo, chain = c(0), component = 0, lag.max = 100, type = "correlation")

par(mfrow = c(1,2))
plot(acfss[[1]][[1]])
plot(acfss[[1]][[2]])
