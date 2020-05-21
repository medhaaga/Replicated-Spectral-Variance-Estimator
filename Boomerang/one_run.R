set.seed(1)
library(cubature)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
source("functions.R")
sourceCpp("lag.cpp")


######### Function calculates CCF with demean available
ccf2 <- function (x, y, lag.max = NULL, type = c("correlation", "covariance"), 
    plot = TRUE, demean = TRUE, na.action = na.fail, ...) 
{
    type <- match.arg(type)
    if (is.matrix(x) || is.matrix(y)) 
        stop("univariate time series only")
    X <- ts.intersect(as.ts(x), as.ts(y))
    colnames(X) <- c(deparse(substitute(x))[1L], deparse(substitute(y))[1L])
    acf.out <- acf(X, lag.max = lag.max, plot = FALSE, type = type, 
        na.action = na.action, demean = FALSE)
    lag <- c(rev(acf.out$lag[-1, 2, 1]), acf.out$lag[, 1, 2])
    y <- c(rev(acf.out$acf[-1, 2, 1]), acf.out$acf[, 1, 2])
    acf.out$acf <- array(y, dim = c(length(y), 1L, 1L))
    acf.out$lag <- array(lag, dim = c(length(y), 1L, 1L))
    acf.out$snames <- paste(acf.out$snames, collapse = " & ")
    if (plot) {
        plot(acf.out, ...)
        return(invisible(acf.out))
    }
    else return(acf.out)
}
#########

A <- 1
B <- 3
C <- 8

samples <- perspective(A, B, C, 1000)
# pdf(file = "contour_plot.pdf")
contour(samples$x, samples$y, samples$z, main="Contour Plot", nlevels = 20, asp = 1)
# filled.contour(samples$x,samples$y,samples$z, color=terrain.colors, main="Contour Plot",)
# dev.off()

start <- matrix(c(C,1,1,C), 2, 2)

nsim <- 1e3
chain1 <- markov.chain(A, B, C, nsim, start[1,])
chain2 <- markov.chain(A, B, C, nsim, start[2,])

par(mfrow = c(1,1))
plot(chain1, xlim = range(c(chain1[,1], chain2[,1])), ylim = range(c(chain1[,2], chain2[,2])), asp = 1)
points(chain2, col = "red")


c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
c2.cen.loc <- scale(chain2, scale = FALSE)

global.mean <- colMeans(rbind(chain1, chain2))

c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
c2.cen <- scale(chain2, center = global.mean, scale =FALSE)
c1.var <- t(c1.cen) %*% (c1.cen)/(nsim)

par(mfrow = c(2, 3))
acf(c1.cen.loc[,1], demean = FALSE, lag.max = 150, ci = 0) #Old ACFS -- look too good to be true
acf(c1.cen.loc[,2], demean = FALSE, lag.max = 150, ci = 0)
ccf(c1.cen.loc[,1], c1.cen.loc[,2], lag.max = 150, ci = 0)


foo1 <- acf(c1.cen[,1], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance") # New Acfs better
foo2 <- acf(c1.cen[,2], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance")
foo3 <- ccf2(c1.cen[,1], c1.cen[,2], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance")
foo1$acf <- foo1$acf/c1.var[1,1]
foo2$acf <- foo2$acf/c1.var[2,2]
foo3$acf <- foo3$acf/sqrt(com.var[1,1] *com.var[2,2])
plot(foo1, ylim = c(0,1))
plot(foo2, ylim = c(0,1))
plot(foo3, ylim = c(-1, 0))




nsim <- 1e4
chain1 <- markov.chain(A, B, C, nsim, start[1,])
chain2 <- markov.chain(A, B, C, nsim, start[2,])

par(mfrow = c(1,1))
plot(chain1, xlim = range(c(chain1[,1], chain2[,1])), ylim = range(c(chain1[,2], chain2[,2])) )
points(chain2, col = "red")

c1.cen.loc <- scale(chain1, scale = FALSE)  ## X_st - bar(X)_s
c2.cen.loc <- scale(chain2, scale = FALSE)

global.mean <- colMeans(rbind(chain1, chain2))

c1.cen <- scale(chain1, center = global.mean, scale =FALSE)
c2.cen <- scale(chain2, center = global.mean, scale =FALSE)
c1.var <- t(c1.cen) %*% (c1.cen)/(nsim)



par(mfrow = c(2, 3))
acf(c1.cen.loc[,1], demean = FALSE, lag.max = 150, ci = 0) #Old ACFS look similar to new acfs now
acf(c1.cen.loc[,2], demean = FALSE, lag.max = 150, ci = 0)
ccf(c1.cen.loc[,1], c1.cen.loc[,2], lag.max = 150, ci = 0)


foo1 <- acf(c1.cen[,1], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance") # New Acfs
foo2 <- acf(c1.cen[,2], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance")
foo3 <- ccf2(c1.cen[,1], c1.cen[,2], demean = FALSE, lag.max = 150, plot = FALSE, type = "covariance")
foo1$acf <- foo1$acf/c1.var[1,1]
foo2$acf <- foo2$acf/c1.var[2,2]
foo3$acf <- foo3$acf/sqrt(com.var[1,1] *com.var[2,2])
plot(foo1, ylim = c(0,1), ylab = "ACF")
plot(foo2, ylim = c(0,1), ylab = "ACF")
plot(foo3, ylim = c(-1, 0), ylab = "ACF")
