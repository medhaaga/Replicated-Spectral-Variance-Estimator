set.seed(130290)
library(mcmcse)
library(mcmc)

source("spatio_temp_function.R")
load("spatio_temp_truth_2e6by1000")
truth <- colMeans(spatio_truth)
data(logit)
level <- .90

## Fitting model with intercept
N <- 1e5

chain <- spatio_temp(N = N, verbose = FALSE)

pdf("spatio_acf_trace_13.pdf")
par(mfrow = c(2,2))
acf(chain[,1], main = expression(paste("Autocorrelation for ", beta[1])))
ccf(chain[,1], chain[,3], main = expression(paste("Cross-correlation for ", beta[1], " and ", beta[3])))
# acf(chain[,1], main = expression(paste("Autocorrelation for ", beta[3])))
plot.ts(chain[,1], ylab = expression(beta[1]),  main = expression(paste("Trace plot for ", beta[1])), xlab = "Iteration")
plot.ts(chain[,3],  ylab = expression(beta[3]), main = expression(paste("Trace plot for ", beta[3])), xlab = "Iteration")
dev.off()
