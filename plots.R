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

#Perspective plot
samples <- perspective(A,C,500)       #from functions.R
persp(samples$x,samples$y,samples$z,col="lightblue",main="Perspective Plot")

#Contour plot
pdf(file = "contour_plot.pdf")
contour(samples$x,samples$y,samples$z,main="Contour Plot")
filled.contour(samples$x,samples$y,samples$z,color=terrain.colors,main="Contour Plot",)
dev.off()

###############################################

A <- 1
C <- 7
freq <- 1e3     #no. of reps for each nsim
start <- matrix(c(2*C,1,1,2*C),2,2)

##################################################################################

#Convergence plots for nsim = 1e3 to nsim=1e5
load(file = "1e3-1e5_RSV_ASV_samp.Rdata")
conv.pts <- seq(1e3,1e6,500)
#Sigma_11
pdf(file = "conv_plot_Sigma_11.pdf")
plot(conv.pts,asv.samp[1,1,], type = "l", col = "red")
lines(conv.pts, rsv.samp[1,1,], col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()

#Sigma_22
pdf(file = "conv_plot_Sigma_22.pdf")
plot(conv.pts,asv.samp[2,2,], type = "l", col = "red")
lines(conv.pts, rsv.samp[2,2,], col="blue")
dev.off()

#Determinant
det.rsv <- rep(0,length(conv.pts))
det.asv <- rep(0,length(conv.pts))
for (i in 1:length(conv.pts)){
  det.rsv[i] <- det(rsv.samp[,,i])
  det.asv[i] <- det(asv.samp[,,i])
}
pdf(file = "conv_plot_det.pdf")
plot(conv.pts,det.asv, type = "l", col="red")
lines(conv.pts,det.rsv, col="blue")
legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
dev.off()

#########################################################################################
#Load the 1000 replications of ASV, RSV and coverage probabilities for each value of nsim
check.pts <- c(1e3,1e4)  #let these be the checkpoints

r <- length(check.pts)

for (i in 1:r){
  #Scatter plot of m Markov chains
  nsim <- check.pts[i]
  chain1 <- markov.chain(A,C,nsim,start[1,])
  chain2 <- markov.chain(A,C,nsim,start[2,])
  #Scatter plot
  pdf(file = paste("scatter_plot_",i,".pdf", sep = ""))
  plot(chain1[,1], chain1[,2], col = "red", xlim = c(-2,2*C), ylim = c(-2,2*C), main = paste("Scatter plot of two Markov chains, n = ", nsim), xlab = "X", ylab = "Y")
  points(chain2[,1], chain2[,2], col = "green")
  dev.off()
  
  #Comparitative density plots of ASV and RSV
  load(file = paste("out",i,".Rdata", sep = ""))
  
  #Sigma_11
  pdf(file = paste("out",i,"_sigma11.pdf", sep = ""))
  plot(density(asv.samp[1,1,]), col="red")
  lines(density(rsv.samp[1,1,]), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()
  
  #Sigma_22
  pdf(file = paste("out",i,"_sigma22.pdf", sep = ""))
  plot(density(asv.samp[2,2,]), col="red")
  lines(density(rsv.samp[2,2,]), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()
  
  #Determinant
  
  det.rsv <- rep(0,freq)
  det.asv <- rep(0,freq)
  for (j in 1:freq){
    det.rsv[j] <- det(rsv.samp[,,j])
    det.asv[j] <- det(asv.samp[,,j])
  }
  
  pdf(file = paste("out",i,"_det.pdf", sep = ""))
  plot(density(det.asv), col="red")
  lines(density(det.rsv), col="blue")
  legend("topright", legend=c("ASV", "RSV"),col=c("red", "blue"), lty=1, cex=1.2)
  dev.off()
}
