####################################
#								   #
######### Simulation 1 #############
#								   #
# MagnoliaNetworkSamplingExample.R #
#								   #
####################################


#Required libraries
require(coda, quiet=TRUE) #mcmc convergence diagnostics
require(mcmcse, quiet=TRUE) #mcmc se calculations
require(network, quiet=TRUE) #base network pacakage
require(sna, quiet=TRUE) #built from network
require(ggplot2, quiet=TRUE) #for graphics
require(plyr, quiet=TRUE) #for data manipulation
require(reshape2, quiet=TRUE) #for data manipulation
require(dplyr, quiet=TRUE) #for data manipulation
require(data.table, quiet=TRUE)  #for data manipulation
require(parallel, quiet=TRUE) #for parallel processing
require(doMC, quiet=TRUE) #for parallel processing
require(xtable, quiet=TRUE) #easy to print tables
require(gridExtra, quiet=TRUE) #for graphics

#####################################################################################################
# Data preperation function
#####################################################################################################

dataPrepHighSchool <-function(data){
  if(class(data)[1]=="network"){
    data.edgelist = as.edgelist.sna(data)
  } else if(class(data)[1]=="matrix"){
    if(which.matrix.type(data)=="adjacency" | which.matrix.type(data)=="edgelist"){
      data.edgelist = as.edgelist.sna(data)
    } else if(which.matrix.type(data)!="adjacency" & which.matrix.type(data)!="edgelist"){
      stop("Error: Matrix must be an adjacency matrix or edgelist matrix")
    }
  } else if (class(data)[1]=="edgelist"){
    data.edgelist = as.edgelist.sna(data)
  } else if(class(data)!="network" & class(data)!="matrix" & data(class)!=c("edgelist","matrix")){
    stop("Error: Improper dataset format. Inputs must be network, adjacency matrix or edgelist")
  }
  
  #True network size
  size = attr(data.edgelist, "n")
  
  #Network size removing isolated users
  nvertices = length(unique(data.edgelist[,1]))
  
  #List of unique nodes
  nodes = sort(unique(data.edgelist[,1]))
  
  #out-degree for each node (same as in-degree for undirected graphs)
  degree.connected = data.frame(node = nodes, degree = rep(0,nvertices))
  
  degree.connected$degree = as.vector(table(data.edgelist[,1]))
  degree.connected$triples = choose(degree.connected$degree, 2)
  degree.connected$triangles = kcycle.census(data.edgelist, mode="graph")$cycle.count[2,-1]
  degree.connected$sex = get.vertex.attribute(data, "Sex")
  degree.connected$grade = get.vertex.attribute(data, "Grade")
  degree.connected$race = get.vertex.attribute(data, "Race")
  
  #create df to use for walks 
  tempdf= data.table(as.data.frame(data.edgelist))
  colnames(tempdf) <- c("snd", "rec", "value")
  tempdf  = setorder(tempdf, snd, rec) 
  connected.df = as.data.frame(tempdf[, .I[.N], by=tempdf$snd])
  connected.df$V0 = c(1, connected.df$V1[-length(connected.df$V1)]+1)
  connected.df = connected.df[,c("tempdf", "V0", "V1")]
  colnames(connected.df) <- c("Node", "Start", "End")
  
  #Gather information for later functions
  return(list(data.edgelist = data.edgelist, nodes = nodes, networksize = nvertices, 
              size = size, degree.connected = degree.connected, 
              connected.df=connected.df, tempdf=tempdf))
}
#####################################################################################################
#####################################################################################################
# Metropolis-Hastings Random Walk Sample Function
#####################################################################################################


metropolisHastingsRandomWalk <- function(nrep, nodes, networksize, data.edgelist, degree.connected, connected.df, tempdf,
                                         updatedmean=TRUE, chains=1, start=NA, cores=1, replace=TRUE){
  
  #Number of steps (nrep) must be greater than or equal to 100
  if(nrep < 100) stop("Number of steps (nrep) must be greater than or equal to 100")
  
  #Starting vector 
  if(sum(is.na(start))==0 & length(start)==chains){
    step1=start
  } else if(sum(is.na(start))!=0 & length(start)==1){
    step1 = sample(nodes, chains, replace=replace)
  } else {
    stop("Starting vector not correct length or has missing values")
  }
  
  if(is.na(start) & replace==TRUE){
    warning("Default settings (start=NA, replace=TRUE) allow starting points to be not unique. For unique points, set replace=FALSE")
  }
  
  #List for each chain, contains one vector of length nrep
  MHRWrandomWalks = rep(list(numeric(nrep)),chains)
  
  #Fill in the first position of each chain with a step1
  for(i in 1:chains){
    MHRWrandomWalks[[i]][1] = step1[i]
  }
  
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  temp = data.frame(data.edgelist)
  MHstep <- function(MHRWrandomWalks){
    k = c(rep(0,nrep))
    prog.i = txtProgressBar(min = 2, max = nrep, style = 3) #progress bar
    for(i in 2:nrep){
      startANDend = as.numeric(connected.df[connected.df$Node==MHRWrandomWalks[[i-1]], c("Start", "End")])
      rowindex = sample(startANDend[1]:startANDend[2], 1)
      MHRWrandomWalks[i] = tempdf$rec[rowindex]
      u = runif(1,0,1)
      k[i-1] = degree.connected$degree[match(MHRWrandomWalks[i-1], degree.connected$node)]
      k[i] = degree.connected$degree[match(MHRWrandomWalks[i], degree.connected$node)]
      
      if ( u <= k[i-1]/k[i] ) {
        MHRWrandomWalks[i] = MHRWrandomWalks[i]
      } else {
        MHRWrandomWalks[i] = MHRWrandomWalks[i-1]
      }
      setTxtProgressBar(prog.i, i)
    }
    # close(prog.i)
    return(MHRWrandomWalks)
  }
  
  MHRWrandomWalkComplete = mclapply(MHRWrandomWalks, MHstep, mc.cores=cores)
  
  MHRWsampleCollection= mclapply(MHRWrandomWalkComplete, function(x){
    degree.connected[match(x, degree.connected$node), -1]}, mc.cores=cores)
  
  return(MHRWsampleCollection)
}



#####################################################################################################
# Gather Data 
#####################################################################################################

library(ergm)
data(faux.magnolia.high)
magnolia = faux.magnolia.high

#Ensure well connected component only 
notwellconnected = which(component.largest(magnolia, connected="weak")==FALSE)
delete.vertices(magnolia, notwellconnected)
dp = dataPrepHighSchool(magnolia)

#get an assortativy measure
gmagnolia = intergraph::asIgraph(magnolia)
igraph::assortativity_degree(gmagnolia)

#Get degree summary of different groups
temp = dp$degree.connected
temp$race = as.factor(temp$race)
temp$race = factor(temp$race, levels(temp$race)[c(6,2,1,3,4,5)])
temp$sex = as.factor(temp$sex)

g.degree.race = ggplot(temp, aes(x=race, y=degree)) +
  geom_boxplot() +
  xlab("race") 

g.degree.grade = ggplot(temp, aes(x=as.factor(grade), y=degree)) +
  geom_boxplot() +
  xlab("grade")
g.degree.sex = ggplot(temp, aes(x=as.factor(sex), y=degree)) +
  geom_boxplot() +
  xlab("sex")

# Plot true distribution
png("DegreeDistByGroup.png")
grid.arrange(g.degree.grade, g.degree.sex, g.degree.race, ncol=3)
dev.off()


#####################################################################################################
# Create sampling shells
#####################################################################################################

# shell metropolis hasting random walk function
magnolia.mhrw.func <- function(nrep, chains, start){
  output = metropolisHastingsRandomWalk(nrep=nrep, 
                                        nodes = dp$nodes, 
                                        networksize = dp$networksize, 
                                        data.edgelist = dp$data.edgelist, 
                                        degree.connected = dp$degree.connected, 
                                        connected.df = dp$connected.df, 
                                        tempdf = dp$tempdf,
                                        updatedmean=FALSE, 
                                        chains=chains, 
                                        start=start, 
                                        replace=TRUE, 
                                        cores=1)
  return(output)
}


# Metropolis Hastings Random Walk estimate function
MHRWestimation <- function(data){
  data$degree = data$degree 
  data$triples = choose(data$degree, 2)
  data$clustercoef = (ifelse(data$triples>0,data$triangles/data$triples,0))
  data$binarysex = ifelse(data$sex=="F", 1, 0)
  data$binaryracewhite = ifelse(data$race=="White", 1, 0)
  return(data)
}


##calculate T^2 statistics for a given sample mean (x), true mean (mean),
##and estimated covariance matrix - sigma
t2.stat <- function(x,mean,sigma,n){
  return(n*t(mean - x) %*% solve(sigma) %*% (mean - x))
}

##FFT function for calculating RSVe on a centered matrix A
mSVEfft <- function (A, b, method = "bartlett")
{
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag2(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))
  FF <- matrix(0, ncol = p, nrow = 2*n)
  FF[1:n,] <- A
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- mvfftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n, ]) / n )
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)
    FF <- FF * matrix(w, nrow = 2*n, ncol = p)
    FF <- fftw_c2r(FF) / (2* n )
    return ((t(A) %*% FF[1:n]) / n )
  }
  
}

confidence_interval <- function(vector, interval) {
  
  vec_sd <- sd(vector)
  n <- length(vector)
  vec_mean <- mean(vector)
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  result <- c("mean" = vec_mean, "lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

