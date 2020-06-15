# source("spatio_temp_function.R")
library(spBayes)
library(parallel)
data("NETemp.dat")
ne.temp <- NETemp.dat

ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]##subset first 2 years (Jan 2000 - Dec. 2002)
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months

coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))

##set starting and priors
p <- 2 #number of regression parameters in each month
starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
"sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
"sigma.eta"=diag(rep(0.01, p)))

tuning <- list("phi"=rep(6, N.t))

priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
"phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
"sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
"tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
"sigma.eta.IW"=list(2, diag(0.001,p)))

##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)
n.samples <- 2e6


beta0 <- 0
beta <- 0
theta <- 0
sigma.eta <- 0
u <- 0

spatio_truth_loop <- function(i)
{
print(i)

chain <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
starting=starting, tuning=tuning, priors=priors, get.fitted =FALSE,
cov.model="exponential", n.samples=n.samples, n.report=1e4, verbose = FALSE)

beta0 <- colMeans(chain$p.beta.0.samples)
beta <- colMeans(chain$p.beta.samples)
theta <- colMeans(chain$p.theta.samples)
sigma.eta <- colMeans(chain$p.sigma.eta.samples[ ,-3])
u <- colMeans(t(chain$p.u.samples))
return(c(beta0, beta, theta, sigma.eta, u))
}

set.seed(130290)
spatio_truth <- mclapply(1:1000, spatio_truth_loop, mc.cores = 5)
spatio_truth <- Reduce(rbind, spatio_truth)
# spatio_truth <- colMeans(spatio_truth)
save(spatio_truth, file = "spatio_temp_truth_2e6by1000")
