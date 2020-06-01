

###################################################
### chunk number 7: mida1
###################################################

library("MCMCpack")
load("mida.rda")

n <- 1e6
m <- 1e3

mida_means <- matrix(0, nrow = m, ncol = 7)
for(i in 1:m)
{
	print(i)
	out <- MCMCpoissonChange(mida~ 1, m=6, c0=13, d0=1,  
		marginal.likelihood=c("Chib95"), mcmc = n, burnin = 0, verbose = 5e5, seed = i)
	mida_means[i, ] <- colMeans(out)
}

save(mida_means, file = "true_mida_1e6_1e3")


# out1 <- MCMCpoissonChange(mida~ 1, m=6, c0=13, d0=1,  
# 		marginal.likelihood=c("Chib95"), mcmc = 1e4, burnin = 0, verbose = 1e6, seed = round(1e3*runif(1)))

# out2 <- MCMCpoissonChange(mida~ 1, m=6, c0=13, d0=1,  
#   marginal.likelihood=c("Chib95"), mcmc = 1e4, burnin = 0, seed = round(1e3*runif(1)), verbose = 1e4)

# par(mfrow = c(2,7))
# for(i in 1:7)
# {
#   plot.ts(out1[,i])
#   plot.ts(out2[,i])
# }

# par(mfrow = c(2,7))
# for(i in 1:7)
# {
#   acf(out1[,i])
#   acf(out2[,i])
# }


# cbind(colMeans(out1), colMeans(out2))

# par(mfrow = c(2,7))
# for(i in 1:7)
# {
#   plot(density(out1[,i]))
# plot(density(out2[,i]))
# }

