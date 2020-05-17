
combined_acf <- function(x, center, chain, component, lag.max = NULL, type = c("covariance", "correlation")){
  chain.cen <- scale(x[[chain]], center = TRUE, scale = FALSE)
  chain.global.cen <- scale(x[[chain]], center = center, scale = FALSE)
  acf.list <- list()
  for (c in 1:length(component))
  {
    j <- component[c]
    normal.acf <- acf(chain.cen[,j], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
    global.acf <- acf(chain.global.cen[,j], type = type, plot = FALSE, demean = FALSE, lag.max = lag.max)
    acf.list[[c]] <- list(normal.acf, global.acf)
  }
  return(acf.list)
}

combined_ccf <- function(x, center, chain, component, lag.max = NULL, type = c("covariance", "correlation")){
  mc.chain <- x[[chain]]
  n <- nrow(mc.chain)
  mu <- colMeans(mc.chain)
  ccf.list <- list()
  i <- component[1]
  j <- component[2]
  mu1 <- mu[i]
  mu2 <- mu[j]
  normal.ccf <- ccf(mc.chain[,i], mc.chain[,j], type = type, plot = FALSE, lag.max = lag.max)
  global.ccf <- ccf(mc.chain[,i], mc.chain[,j], type = type, plot = FALSE, lag.max = lag.max)
  for (k in 1:150){
    res1 <- (n - k)*(mu1 - center[i])*(mu2 - center[j])
    res2 <- (n - k)*(mu1 - center[i])*(mu2 - center[j])
    for (t in 1:k){
      res1 = res1 + (mu1 - center[i])*(mu2 - mc.chain[t,j]) + (mu2 - center[j])*(mu1 - mc.chain[n-k+t,i])
      res2 = res2 + (mu1 - center[i])*(mu2 - mc.chain[n-k+t,j]) + (mu2 - center[j])*(mu1 - mc.chain[t,i])
    }
    global.ccf$acf[k] = global.ccf$acf[k] + res1/n
    global.ccf$acf[301-k+1] = global.ccf$acf[301-k+1] +res2/n
  }
  global.ccf$acf[151] = global.ccf$acf[151] + (mu1 - center[i])*(mu2 - center[j])
  ccf.list <- list(normal.ccf, global.ccf)
  return(ccf.list)
}

