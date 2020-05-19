set.seed(1)
library(fields)
library(scatterplot3d)
library(cubature)
library(graphics)
library(pracma)
library(Rcpp)
library(RcppArmadillo)
library(fftwtools)
library(mcmcse)
library(stats)
sourceCpp("lag.cpp")
source("functions.R")


###############################################
#####Visualization of distribution############
###############################################
# Uncomment to visualize through contour plots
#Contour plot
A <- 2
B <- 9
C <- 7
samples <- perspective(A, B, C, 500)
pdf(file = paste("Out/", A, B, C, "/3d_density_plot.pdf", sep = "_"), height = 4)
par(mfrow = c(1,2))
contour(samples$x,samples$y,samples$z,main="Contour Plot")
#filled.contour(samples$x,samples$y,samples$z,color=terrain.colors,main="Contour Plot",)
persp(samples$x,samples$y,samples$z,main="Perspective Plot", col = "springgreen", shade = .5, theta = -45, phi = 0, xlab = "x", ylab = "y", zlab = "Density")
dev.off()

###############################################

m = 5
check.pts <- c(1e3, 2e3, 5e3, 1e4, 2e4)

r <- length(check.pts)
freq <- 1e3
c.prob <- .95
min <- 5e2
max <- 5e4
conv.pts <- seq(min, max, 100)

  A <- 2
  B <- 9
  C <- 7
  start <- matrix(0, nrow = m, ncol = 2)  #only depends on C

  for(i in 1:floor(m/2)){
    start[i,] <- c(0, C*(2^(2-i)))
    start[m-i+1,] <- c(C*(2^(2-i)), 0)
  }

  ## 1.) Running plots
  load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

  #### 1a.) Sigma_11
  pdf(file = paste(paste("Out/", A, B, C, "/run_plots", m, sep = "_"), ".pdf", sep = ""), height = 3)
  par(mfrow = c(1,3))

  plot(conv.pts,asv.samp[1,1,], type = "l", col = "red", main = paste("Variance of x-component"), xlab = "Simulation size", ylab = "Variance")
  lines(conv.pts, rsv.samp[1,1,], col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

  plot(conv.pts,asv.samp[2,2,], type = "l", col = "red", main = "Variance of y-component", xlab = "Simulation size", ylab = "Variance")
  lines(conv.pts, rsv.samp[2,2,], col="blue")

  #### 1c.) Determinant
  det.rsv <- rep(0,length(conv.pts))
  det.asv <- rep(0,length(conv.pts))

  for (i in 1:length(conv.pts)){
    det.rsv[i] <- det(rsv.samp[,,i])
    det.asv[i] <- det(asv.samp[,,i])
  }

  plot(conv.pts,det.asv, type = "l", col="red", main = paste("Determinant"), xlab = "Simulation size", ylab = "Determinant")
  lines(conv.pts,det.rsv, col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()

  #### 1d.) ESS

  pdf(file = paste("Out/", A, B, C, "/run_plot_ess_2_5.pdf", sep = "_"), height = 4)

  m<-2
  load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))

  par(mfrow = c(1,2))
  mean.asv <- rep(0, 96)
  mean.rsv <- rep(0, 96)
  #lower.asv <- rep(0, length(conv.pts))
  #lower.rsv <- rep(0, length(conv.pts))
  #upper.asv <- rep(0, length(conv.pts))
  #upper.rsv <- rep(0, length(conv.pts))
  for (i in 1:96){
    asv <- confidence_interval(ess.asv.samp[[i]], .9)
    rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
    mean.asv[i] <- asv[1]
    mean.rsv[i] <- rsv[1]
    #lower.asv[i] <- asv[2]
    #lower.rsv[i] <- rsv[2]
    #upper.asv[i] <- asv[3]
    #upper.rsv[i] <- rsv[3]
  }
  plot(conv.pts[1:96], mean.asv, type = "l", col = "red", main = "Two parallel chains", xlab = "Simulation size", ylab = "ESS/mn")
  lines(conv.pts[1:96], mean.rsv, type = "l", col = "blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  m = 5
  load(file = paste(paste("Out/", A, B, C, "/conv_data", m, min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))
  for (i in 1:96){
     asv <- confidence_interval(ess.asv.samp[[i]], .9)
     rsv <- confidence_interval(ess.rsv.samp[[i]], .9)
     mean.asv[i] <- asv[1]
     mean.rsv[i] <- rsv[1]
     #lower.asv[i] <- asv[2]
     #lower.rsv[i] <- rsv[2]
     #upper.asv[i] <- asv[3]
     #upper.rsv[i] <- rsv[3]
   }
  plot(conv.pts[1:96], mean.asv, col = "red", type = "l", main = "Five parallel chains", xlab = "Simulation size", ylab = "ESS/mn")
  lines(conv.pts[1:96], mean.rsv, col = "blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()


  ## 2.) Density plots for Sigma11, Sigma22 and determinant.
  ##     Loop over check points for different values of nsim

  for (j in 1:r){

    #### 2a.) Scatter plot of m Markov chains
    nsim <- check.pts[j]
    chain <- array(0, dim = c(nsim,2,m))

    for (p in 1:m){
      chain[,,p] <- markov.chain(A, B, C, nsim, start[p,])
    }

    pdf(file = paste(paste("Out/", A, B, C, "/scatter_plot", m, nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    plot(chain[,1,1], chain[,2,1], col = 2, xlim = c(-2,2*C), ylim = c(-2,2*C), main = paste("Scatter plot of m = ", m,  "Markov chains, n = ", nsim, ", A = ", A, ", B = ", B, ", C = ", C), xlab = "X", ylab = "Y")
    for (p in 2:m){
      points(chain[,1,p], chain[,2,p], col = p+1)
    }

    dev.off()

    pdf(file = paste(paste("Out/", A, B, C, "/trace", m, nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    par(mfrow = c(m,2))
    for(p in 1:m){
      plot.ts(chain[,1,p], main = paste("x component of chain -", p))
      plot.ts(chain[,2,p], main = paste("y component of chain -", p))
    }
    dev.off()

    #Comparitative density plots of ASV and RSV
    load(file = paste(paste("Out/", A, B, C, "/out", m, nsim, A, B, C, sep = "_"), ".Rdata", sep = ""))

    #### 2b.) Sigma_11
    pdf(file = paste(paste("Out/", A, B, C, "/densities",m, nsim, sep = "_"), ".pdf", sep = ""), height = 3)
    par(mfrow = c(1,3))
    plot(density(asv.samp[1,1,]), col="red", main = "Variance of x-component", xlab = "Variance")
    lines(density(rsv.samp[1,1,]), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

    #### 2c.) Sigma_22

    plot(density(asv.samp[2,2,]), col="red", main = "Variance of y-component", xlab = "Variance")
    lines(density(rsv.samp[2,2,]), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)

    #### 2d.) Determinant

    det.rsv <- rep(0,freq)
    det.asv <- rep(0,freq)
    for (k in 1:freq){
      det.rsv[k] <- det(rsv.samp[,,k])
      det.asv[k] <- det(asv.samp[,,k])
    }

    plot(density(det.asv), col="red", main = "Determinant", xlab = "Determinant")
    lines(density(det.rsv), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
    dev.off()

    print(paste("Coverage probabilities for nsim =  ", nsim, ", m = ", m, ", A = ", A, ", B = ", B, ", C = ", C, " are: "))
    print(paste("ASV : ", mean(asv.coverage), "RSV: ", mean(rsv.coverage)))
  }
#}

