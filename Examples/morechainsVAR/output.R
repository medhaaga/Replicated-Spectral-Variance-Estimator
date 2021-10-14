
library(multichainACF)
load(file = "Out/var-100_chains.Rdata")
#######################################################
############### ACF and G-ACF #########################
#######################################################

component <- 1
lag.max <- 40

m <- 101
lag1 <- 1
lag2 <- lag.max

lacf_twolags <- matrix(0, nrow = 2, ncol = m)
gacf_twolags <- matrix(0, nrow = 2, ncol = m)

nsims <- c(1e3, 1e4)

for (n in 1:2)
{
  ncrop <- nsims[n]
  
  for (t in 2:m)
  {
    print(t)
    nchains <- t
    
    
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
  
  pdf(file = paste("Out/ACFvsm_lag1_", ncrop, ".pdf", sep = ""), height=5, width=5)
  plot(2:m, lacf_twolags[1,-1], type = "l", ylim = c(.97,1), col = "blue", ylab = "ÄCF", xlab = "m", main = paste("n =", ncrop))
  lines(2:m, gacf_twolags[1,-1], col = "orange")
  abline(h = true.acf[1,1,lag1+lag.max+1]/true.acf[1,1,lag.max+1], col = "red", lwd=1)
  legend("bottomright", legend = c("A-ACF", "G-ACF", "Truth"), col = rep(c("blue", "orange", "red"), 1), lty=c(1,1,1), cex=.7)
  dev.off()


  pdf(file = paste("Out/ACFvsm_lag40_", ncrop, ".pdf", sep = ""), height=5, width=5)
  plot(2:m, lacf_twolags[2,-1], type = "l", ylim = c(.6,.99), col = "blue", ylab = "ÄCF", xlab = "m", main = paste("n =", ncrop))
  lines(2:m, gacf_twolags[2,-1], col = "orange", lty=1)
  abline(h = true.acf[1,1,lag2+lag.max+1]/true.acf[1,1,lag.max+1], col = "red", lwd=1, lty=1)
  legend("bottomright", legend = c("A-ACF", "G-ACF", "Truth"), col = rep(c("blue", "orange", "red"), 1), lty=c(1,1,1), cex=.7)
  dev.off() 
}

