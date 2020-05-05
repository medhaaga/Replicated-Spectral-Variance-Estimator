source("functions.R")
sourceCpp("lag.cpp")

##creates freq=1000 replications of ASv and RSV for each value of nsim from check.pts

create.output <- function(A,C,m,check.pts,freq,start,c.prob,T.mean){
  critical <- qchisq(c.prob, df=2)
  for (i in 1:length(check.pts)){
    nsim <- check.pts[i]
    asv.samp <- array(0, dim = c(2,2,freq))
    rsv.samp <- array(0, dim = c(2,2,freq))
    asv.coverage <- rep(0,freq)
    rsv.coverage <- rep(0,freq)
    
    for (j in 1:freq){
      chain <- array(0,dim = c(nsim,2,m))
      sve <- array(0, dim = c(2,2,m))
      rsve <- array(0, dim = c(2,2,m))
      b <- rep(0,m)
      for (k in 1:m){
        chain[,,k] <- markov.chain(A,C,nsim,start[k,])
        b[k] <- batchSize(chain[,,k], method = "bartlett")
      }
      b.avg <- mean(b)
      global.mean <- apply(chain,2,mean)
      for (k in 1:m){
        chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
        sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
        chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
        rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
      }
      asv.samp[,,j] <- apply(sve, c(1,2), mean)
      rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
      
      chi.sq.asv <- ((nsim-2)/(nsim-1))*t2.stat(global.mean,T.mean,asv.samp[,,j],nsim)
      chi.sq.rsv <- ((nsim-2)/(nsim-1))*t2.stat(global.mean,T.mean,rsv.samp[,,j],nsim)
      if (chi.sq.asv <= critical) {asv.coverage[j]=1}
      if (chi.sq.rsv <= critical) {rsv.coverage[j]=1}
    }
    save(asv.samp,rsv.samp,asv.coverage,rsv.coverage, file = paste("out",i,".Rdata", sep = ""))
  }
  
}

######################################################################
A=1
C=7
m=2
start <- matrix(c(2*C,1,1,2*C), 2, 2)
T.mean <- true.mean(A,C)
#sims for plotting densities and calculating coverage

check.pts <- c(1e3, 1e4)
freq <- 1e3
c.prob <- .95
create.output(A,C,m,check.pts,freq,start,c.prob,T.mean)

#######################################################################
#for convergence plots
min <- 1e3
max <- 1e6
master.chain <- array(0,dim = c(max,2,m))
#create m markov chains of max iterations
for (k in 1:m){
  master.chain[,,k] <- markov.chain(A,C,max,start[k,])
}
conv.pts <- seq(min, max,500)

l <- length(conv.pts)
asv.samp <- array(0, dim = c(2,2,l))
rsv.samp <- array(0, dim = c(2,2,l))
for (j in 1:l){
  nsim = conv.pts[j]
  chain <- master.chain[1:nsim,,]
  sve <- array(0, dim = c(2,2,m))
  rsve <- array(0, dim = c(2,2,m))
  b <- rep(0,m)
  for (k in 1:m){
    b[k] <- batchSize(chain[,,k], method = "bartlett")
  }
  b.avg <- mean(b)
  global.mean <- apply(chain,2,mean)
  for (k in 1:m){
    chain.cen.loc <- scale(chain[,,k], scale = FALSE)  ## X_st - bar(X)_s
    sve[,,k] <- mSVEfft(A = chain.cen.loc, b = b.avg)
    chain.cen <- scale(chain[,,k], center = global.mean, scale =FALSE)
    rsve[,,k] <- mSVEfft(A = chain.cen, b = b.avg)
  }
  asv.samp[,,j] <- apply(sve, c(1,2), mean)
  rsv.samp[,,j] <- apply(rsve, c(1,2), mean)
  print(paste("Done with j=", j))
}
save(asv.samp,rsv.samp, file = "1e3-1e5_RSV_ASV_samp.Rdata")

####################################################################33
