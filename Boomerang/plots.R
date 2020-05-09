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
#A <- 2
#B <- 10
#C <- 8
#samples <- perspective(A, B, C, 500)
#pdf(file = "contour_plot.pdf")
#contour(samples$x,samples$y,samples$z,main="Contour Plot")
#filled.contour(samples$x,samples$y,samples$z,color=terrain.colors,main="Contour Plot",)
#dev.off()

###############################################

params <- matrix(c(1,2,7,1,8,9,1,9,9,2,6,7,2,8,7,2,9,7,2,10,8), nrow = 7, ncol=3, byrow = TRUE)
m = 2
check.pts <- c(1e3, 5e3, 1e4, 5e4, 7e4, 1e5)
conv.pts <- seq(1e3, 1e5, 500)
r <- length(check.pts)
freq <- 1e3
c.prob <- .95
min <- 1e3
max <- 1e5

#loop over all parameters A, B, C
for (t in 1:7){
  A <- params[t,1]
  B <- params[t,2]
  C <- params[t,3]
  start <- matrix(c(2*C,1,1,2*C), 2, 2)
  print(paste("For parameters A, B, C = ", A, B, C, " respectively."))
  
  ## 1.)Covergence plots
  load(file = paste(paste("conv_data", min, max, A, B, C, sep = "_"), ".Rdata", sep = ""))
  
  #### 1a.) Sigma_11
  pdf(file = paste(paste("conv_plot_Sig11", A, B, C, sep = "_"), ".pdf", sep = ""))
  plot(conv.pts,asv.samp[1,1,], type = "l", col = "red", main = paste("Convergence plot of Sigma_11, A = ", A, ", B = ", B, ", C = ", C ), xlab = "Simulation size")
  lines(conv.pts, rsv.samp[1,1,], col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()
  
  #### 1b.) Sigma_22
  pdf(file = paste(paste("conv_plot_Sig22", A, B, C, sep = "_"), ".pdf", sep = ""))
  plot(conv.pts,asv.samp[2,2,], type = "l", col = "red", main = paste("Convergence plot of Sigma_22, A = ", A, ", B = ", B, ", C = ", C ), xlab = "Simulation size")
  lines(conv.pts, rsv.samp[2,2,], col="blue")
  dev.off()
  
  #### 1c.) Determinant
  det.rsv <- rep(0,length(conv.pts))
  det.asv <- rep(0,length(conv.pts))
  
  for (i in 1:length(conv.pts)){
    det.rsv[i] <- det(rsv.samp[,,i])
    det.asv[i] <- det(asv.samp[,,i])
  }
  
  pdf(file = paste(paste("conv_plot_det", A, B, C, sep = "_"), ".pdf", sep = ""))
  plot(conv.pts,det.asv, type = "l", col="red", main = paste("Convergence plot of determinant, A = ", A, ", B = ", B, ", C = ", C ), xlab = "Simulation size")
  lines(conv.pts,det.rsv, col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()
  ## 2.) Density plots for Sigma11, Sigma22 and determinant.
  ##     Loop over check points for different values of nsim
  
  for (j in 1:r){
    
    #### 2a.) Scatter plot of m Markov chains
    nsim <- check.pts[j]
    chain1 <- markov.chain(A, B, C, nsim, start[1,])
    chain2 <- markov.chain(A, B, C, nsim, start[2,])
    
    
    pdf(file = paste(paste("scatter_plot",nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    plot(chain1[,1], chain1[,2], col = "red", xlim = c(-2,2*C), ylim = c(-2,2*C), main = paste("Scatter plot of two Markov chains, n = ", nsim, ", A = ", A, ", B = ", B, ", C = ", C), xlab = "X", ylab = "Y")
    points(chain2[,1], chain2[,2], col = "green")
    dev.off()
    
    #Comparitative density plots of ASV and RSV
    load(file = paste(paste("out",nsim,A,B,C, sep = "_"),".Rdata", sep = ""))
    
    #### 2b.) Sigma_11
    pdf(file = paste(paste("density_S11",nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    plot(density(asv.samp[1,1,]), col="red", main = paste("Density plot for Sigma_11, nsim = ", nsim, ", A = ", A, ", B = ", B, ", C = ", C))
    lines(density(rsv.samp[1,1,]), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
    dev.off()
    
    #### 2c.) Sigma_22
    pdf(file = paste(paste("density_S22",nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    plot(density(asv.samp[2,2,]), col="red", main = paste("Density plot for Sigma_22, nsim = ", nsim, ", A = ", A, ", B = ", B, ", C = ", C))
    lines(density(rsv.samp[2,2,]), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
    dev.off()
    
    #### 2d.) Determinant
    
    det.rsv <- rep(0,freq)
    det.asv <- rep(0,freq)
    for (j in 1:freq){
      det.rsv[j] <- det(rsv.samp[,,j])
      det.asv[j] <- det(asv.samp[,,j])
    }
    
    pdf(file = paste(paste("density_det",nsim, A, B, C, sep = "_"), ".pdf", sep = ""))
    plot(density(det.asv), col="red", main = paste("Density plot for determinant, nsim = ", nsim, ", A = ", A, ", B = ", B, ", C = ", C))
    lines(density(det.rsv), col="blue")
    legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
    dev.off()
  }
}

