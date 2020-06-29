
############### RAM transition used at each iteration
############### target function is defined as a function "target"
############### See line 156 for example.

ram.transition <- function(current.x, current.z, var.cov, epsilon) {

  x.c <- current.x
  z.c <- current.z
  x.dim <- length(x.c)
  x.c.den <- target(x.c)
  z.c.den <- target(z.c)
  accept <- 0

  # downhill
  x.p1 <- mvrnorm(1, mu = x.c, Sigma = var.cov)
  x.p1.den <- target(x.p1)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + epsilon) - log(x.p1.den + epsilon)) {
    x.p1 <- mvrnorm(1, mu = x.c, Sigma = var.cov)
    x.p1.den <- target(x.p1)
    N.d <- N.d + 1
  }

  # uphill
  x.p2 <- mvrnorm(1, mu = x.p1, Sigma = var.cov)
  x.p2.den <- target(x.p2)
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + epsilon) - log(x.p1.den + epsilon)) {
    x.p2 <- mvrnorm(1, mu = x.p1, Sigma = var.cov)
    x.p2.den <- target(x.p2)
    N.u <- N.u + 1
  }

  # downhill for auxiliary variable z
  z <- mvrnorm(1, mu = x.p2, Sigma = var.cov)
  z.den <- target(z)
  N.z <- 1  
  while (-rexp(1) > log(x.p2.den + epsilon) - log(z.den + epsilon)) {
    z <- mvrnorm(1, mu = x.p2, Sigma = var.cov)
    z.den <- target(z)
    N.z <- N.z + 1  
  }


  # accept or reject the proposal
  min.nu <- min(1, (x.c.den + epsilon) / (z.c.den + epsilon))
  min.de <- min(1, (x.p2.den + epsilon) / (z.den + epsilon))
  l.mh <- log(x.p2.den) - log(x.c.den) + log(min.nu) - log(min.de)

  if (l.mh > -rexp(1)) {
    x.c <- x.p2
    z.c <- z
    accept <- 1
  }

  c(x.c, z.c, accept, N.d, N.u, N.z)

}

############### RAM algorithm

ram <- function(x.initial, z.initial, var.cov, sample.size, burn.size, 
                epsilon = 10^(-308), one.time.adaptation = FALSE) {
  
  print(Sys.time())
  n.total <- sample.size + burn.size
  x.dim <- length(x.initial)
  accept <- rep(0, n.total)

  N.d <- rep(NA, n.total)
  N.u <- rep(NA, n.total)
  N.z <- rep(NA, n.total)

  if (x.dim == 1) {
    out <- rep(NA, n.total)
  } else {
    out <- matrix(NA, nrow = n.total, ncol = x.dim)
  }

  x.t <- x.initial
  z.t <- z.initial
  var.t <- var.cov

  for (i in 1 : n.total) {

    temp <- ram.transition(x.t, z.t, var.t, epsilon)

    if (x.dim == 1) {
      x.t <- out[i] <- temp[1]
      z.t <- temp[2]
      accept[i] <- temp[3]
      N.d[i] <- temp[4]
      N.u[i] <- temp[5]
      N.z[i] <- temp[6]
    } else {
      x.t <- out[i, ] <- temp[1 : x.dim]
      z.t <- temp[(x.dim + 1) : (2 * x.dim)]
      accept[i] <- temp[2 * x.dim + 1]
      N.d[i] <- temp[2 * x.dim + 2]
      N.u[i] <- temp[2 * x.dim + 3]
      N.z[i] <- temp[2 * x.dim + 4]
    }   

    if (one.time.adaptation == TRUE) {
      if (i == burn.size) {
        if (x.dim == 1) {
          var.t <- var(out[1 : i])
        } else {
          var.t <- var(out[1 : i, ])
        }
      }
    }

  }

  print(Sys.time())

  if (x.dim == 1) {
    out <- list(x = out[-c(1 : burn.size)],
                accept = accept[-c(1 : burn.size)],
                N.d = N.d,
                N.u = N.u,
                N.z = N.z)
  } else {
    out <- list(x = out[-c(1 : burn.size), ],
                accept = accept[-c(1 : burn.size)],
                N.d = N.d,
                N.u = N.u,
                N.z = N.z)
  }

  out  

}

############################ Example 1: A mixture of 20 bivariate Gaussian distributions

mean.mat <- matrix(c(2.18, 5.76, 8.67, 9.59, 4.24, 8.48, 8.41, 1.68, 3.93, 8.82, 
                     3.25, 3.47, 1.70, 0.50, 4.59, 5.60, 6.91, 5.81, 6.87, 5.40, 
                     5.41, 2.65, 2.70, 7.88, 4.98, 3.70, 1.14, 2.39, 8.33, 9.50, 
                     4.93, 1.50, 1.83, 0.09, 2.26, 0.31, 5.54, 6.86, 1.69, 8.11), 
                   ncol = 2, byrow = TRUE)

# If you have not installed the packages below,
# install.packages("mnormt")
library(mnormt)
# install.packages("MASS")
library(MASS)

############### Case (a)

sigma <- 0.1
weight <- 0.05
varCov <- sigma^2 * diag(2)
target <- function(a) {
  sum(sapply(1 : dim(mean.mat)[1], 
             function(t) { weight * dmnorm(a, mean = mean.mat[t, ], varcov = varCov) }))
}

res.comb <- vector("list", 20)
sc <- 4 * diag(2)

for (j in 1 : 20) {
  res.comb[[j]] <- ram(x.initial = runif(2), z.initial = runif(2), var.cov = sc,
                       sample.size = 100, burn.size = 100)
  # The sample size used in the article is "sample.size = 50000" and "burn.size = 25000".
}

e_x1 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x1[i] <- mean(res.comb[[i]]$x[, 1])
}
mean(e_x1); sd(e_x1)

e_x2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x2[i] <- mean(res.comb[[i]]$x[, 2])
}
mean(e_x2); sd(e_x2)

e_x1_2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x1_2[i] <- mean(res.comb[[i]]$x[, 1]^2)
}
mean(e_x1_2); sd(e_x1_2)

e_x2_2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x2_2[i] <- mean(res.comb[[i]]$x[, 2]^2)
}
mean(e_x2_2); sd(e_x2_2)

######## The number of target density evaluations

N.u.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.u[-c(1 : 25000)]) })
N.d.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.d[-c(1 : 25000)]) })
N.z.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.z[-c(1 : 25000)]) })
mean(N.u.summary) + mean(N.d.summary) + mean(N.z.summary)

############### Case (b)

tau2 <- sapply(1 : 20, function(j) { sqrt(sum((mean.mat[j, ] - c(5, 5))^2)) / 20 })
weight <- 1 / sapply(1 : 20, function(j) { sqrt(sum((mean.mat[j, ] - c(5, 5))^2)) })
weight <- weight / sum(weight)
target <- function(a) {
  sum(sapply(1 : dim(mean.mat)[1], function(t) {
             weight[t] * dmnorm(a, mean = mean.mat[t, ], varcov = tau2[t] * diag(2))
            }))
}

res.comb <- vector("list", 20)
sc <- 3.5 * diag(2)

for (j in 1 : 20) {
  res.comb[[j]] <- ram(x.initial = runif(2), z.initial = runif(2), var.cov = sc,
                       sample.size = 100, burn.size = 100)
  # The sample size used in the article is "sample.size = 50000" and "burn.size = 25000".
}

e_x1 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x1[i] <- mean(res.comb[[i]]$x[, 1])
}
mean(e_x1); sd(e_x1)

e_x2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x2[i] <- mean(res.comb[[i]]$x[, 2])
}
mean(e_x2); sd(e_x2)

e_x1_2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x1_2[i] <- mean(res.comb[[i]]$x[, 1]^2)
}
mean(e_x1_2); sd(e_x1_2)

e_x2_2 <- rep(NA, 20)
for (i in 1 : 20) {
  e_x2_2[i] <- mean(res.comb[[i]]$x[, 2]^2)
}
mean(e_x2_2); sd(e_x2_2)

######## The number of target density evaluations

N.u.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.u[-c(1 : 25000)]) })
N.d.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.d[-c(1 : 25000)]) })
N.z.summary <- sapply(1 : 20, function(j) { mean(res.comb[[j]]$N.z[-c(1 : 25000)]) })
mean(N.u.summary) + mean(N.d.summary) + mean(N.z.summary)

############################ Example 2: High dimensional behavior

############### Metropolis

Metro <- function(initial.x = 8, var.cov, n = 100, burn = 100, 
                  one.time.adaptation = FALSE) {

  print(Sys.time())

  n.total <- n + burn
  x.dim <- length(initial.x)
  if (x.dim == 1) {
    out <- rep(NA, n.total)
  } else {
    out <- matrix(NA, nrow = n.total, ncol = x.dim)
  }
  x.t <- initial.x
  var.t <- var.cov
  accept <- rep(0, n.total)

  for (i in 1 : n.total) {

    x.p <- rmnorm(1, mean = x.t, varcov = var.t)
    l.metro <- log(target(x.p)) - log(target(x.t))

    if (l.metro >  -rexp(1)) {
      x.t <- x.p
      accept[i] <- 1
    }

    if (x.dim == 1) {
      out[i] <- x.t
    } else {
      out[i, ] <- x.t
    }

    if (one.time.adaptation == TRUE) {
      if (i == burn) {
        if (x.dim == 1) {
          var.t <- var(out[1 : i])
        } else {
          var.t <- var(out[1 : i, ])
        }
      }
    }
  }

  print(Sys.time())

  if (x.dim == 1) {
    out <- list(x = out[-c(1 : burn)],
                accept = accept[-c(1 : burn)])
  } else {
    out <- list(x = out[-c(1 : burn), ],
                accept = accept[-c(1 : burn)])
  }

  out

}

############### Parallel tempering

pt <- function(initial.x = 8, temp, varcov, swap.rate = 1, 
               n = 10, burn = 10, one.time.adaptation = FALSE) {

  print(Sys.time())

  n.total <- n + burn
  x.dim <- length(initial.x)
  n.temp <- length(temp)
  out <- matrix(NA, nrow = n.total, ncol = x.dim)
  den.caching <- rep(1, n.temp)

  x.t <- rep(initial.x, n.temp)
  accept <- matrix(0, nrow = n.total, ncol = n.temp)
  swap.accept <- rep(0, n.total)

  for (i in 1 : n.total) {
    x.p <- x.t + as.vector(t(rmnorm(n.temp, mean = rep(0, x.dim), varcov = varcov)))
    x.p.m <- matrix(x.p, ncol = n.temp)
    x.t.m <- matrix(x.t, ncol = n.temp)
    x.p.den <- sapply(1 : n.temp, function(j) {
      pt.target(temp[j], x.p.m[, j])
    })
    x.t.den <- sapply(1 : n.temp, function(j) {
      pt.target(temp[j], x.t.m[, j])
    })
    l.metro <- x.p.den - x.t.den
    den.caching <- x.t.den
    accept.reject <- l.metro > -rexp(n.temp)
    x.t.m[, accept.reject] <- x.p.m[, accept.reject]
    den.caching[accept.reject] <- x.p.den[accept.reject]
    accept[i, ] <- accept.reject
 
    if (rbinom(1, 1, swap.rate)) {
      swap.temp <- sample(1 : (n.temp - 1), 1)
      swap <- c(swap.temp, swap.temp + 1)
      l.mh <- pt.target(temp[swap[1]], x.t.m[, swap[2]]) + 
              pt.target(temp[swap[2]], x.t.m[, swap[1]]) -
              den.caching[swap[1]] - den.caching[swap[2]]

      if (l.mh > -rexp(1)) {
        temp.swap <- x.t.m[, swap[1]]
        x.t.m[, swap[1]] <- x.t.m[, swap[2]]
        x.t.m[, swap[2]] <- temp.swap
        swap.accept[i] <- 1
      }

    }

    x.t <- as.numeric(x.t.m)
    out[i, ] <- x.t[1 : x.dim]

    if (one.time.adaptation == TRUE) {
      if (i == burn) {
        if (x.dim == 1) {
          var.t <- var(out[1 : i])
        } else {
          var.t <- var(out[1 : i, ])
        }
      }
    }
  }

  print(Sys.time())
  out <- list(x = out[-c(1 : burn), ],
              swap.accept = swap.accept[-c(1 : burn)],
              accept = accept[-c(1 : burn), ])
  out  

}

############### Target distribution

mu1 <- c(10, 10, 10, rep(c(0, 10), 4))
mu2 <- c(0, 0, 0, rep(c(10, 0), 4))
mu3 <- c(10, 0, 10, rep(c(0, 10), 4))
mu4 <- c(0, 10, 10, rep(c(0, 10), 4))
mu5 <- c(0, 0, 10, rep(c(0, 10), 4))
mu6 <- c(0, 10, 0, rep(c(10, 0), 4))
mu7 <- c(10, 0, 0, rep(c(10, 0), 4))
mu8 <- c(10, 10, 0, rep(c(10, 0), 4))

target <- function(a) {
  exp(sum(dnorm(a, mean = mu1[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu2[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu3[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu4[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu5[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu6[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu7[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8 +
  exp(sum(dnorm(a, mean = mu8[1 : test.dim], sd = rep(1, test.dim), log = TRUE))) / 8
}

############### Dimension d

test.dim <- 3
# Please set "test.dim" to 5, 7, 9, or 11 for the results in the other dimension.

######## Setting the initial jumping covariance matrix

td <- c(0, 10)
try <- vector("list", 2)

for (i in 1 : 2) {
  try[[i]] <- Metro(initial.x = rep(td[i], test.dim), 
                    var.cov = 2.38^2 / test.dim * diag(test.dim),
                    n = 5000, burn = 0)
}

x.save <- NULL

for (i in 1 : 2) {
  x.save <- rbind(x.save, try[[i]]$x)
}
var.c <- cov(x.save)

######## RAM

ram.comb <- vector("list", 10)

for (j in 1 : 10) {
  ram.comb[[j]] <- ram(rep(0, test.dim), z.initial = rep(0, test.dim), var.cov = var.c, 
                       sample.size = 100, burn.size = 100, one.time.adaptation = TRUE)
}
# The sample size used in the article for all the dimensions is 
# "sample.size = 300000" and "burn.size = 200000",

# Summary
ram.comb.accept <- rep(NA, 10)
d000 <- rep(NA, 10); d100 <- rep(NA, 10); d010 <- rep(NA, 10)
d001 <- rep(NA, 10); d110 <- rep(NA, 10); d101 <- rep(NA, 10)
d011 <- rep(NA, 10); d111 <- rep(NA, 10); avg.err <- rep(NA, 10)
N.d <- rep(NA, 10); N.u <- rep(NA, 10); N.z <- rep(NA, 10)

for (j in 1 : 10) {
  ram.comb.accept[j] <- mean(ram.comb[[j]]$accept)
  d000[j] <- sum(ram.comb[[j]]$x[, 1] < 5 & ram.comb[[j]]$x[, 2] < 5 & ram.comb[[j]]$x[, 3] < 5)
  d100[j] <- sum(ram.comb[[j]]$x[, 1] > 5 & ram.comb[[j]]$x[, 2] < 5 & ram.comb[[j]]$x[, 3] < 5)
  d010[j] <- sum(ram.comb[[j]]$x[, 1] < 5 & ram.comb[[j]]$x[, 2] > 5 & ram.comb[[j]]$x[, 3] < 5)
  d001[j] <- sum(ram.comb[[j]]$x[, 1] < 5 & ram.comb[[j]]$x[, 2] < 5 & ram.comb[[j]]$x[, 3] > 5)
  d110[j] <- sum(ram.comb[[j]]$x[, 1] > 5 & ram.comb[[j]]$x[, 2] > 5 & ram.comb[[j]]$x[, 3] < 5)
  d101[j] <- sum(ram.comb[[j]]$x[, 1] > 5 & ram.comb[[j]]$x[, 2] < 5 & ram.comb[[j]]$x[, 3] > 5)
  d011[j] <- sum(ram.comb[[j]]$x[, 1] < 5 & ram.comb[[j]]$x[, 2] > 5 & ram.comb[[j]]$x[, 3] > 5)
  d111[j] <- sum(ram.comb[[j]]$x[, 1] > 5 & ram.comb[[j]]$x[, 2] > 5 & ram.comb[[j]]$x[, 3] > 5)
  avg.err[j] <- mean(abs(d000[j] / 3e5 - 1 / 8) + abs(d100[j] / 3e5 - 1 / 8) + 
                     abs(d010[j] / 3e5 - 1 / 8) + abs(d001[j] / 3e5 - 1 / 8) + 
                     abs(d110[j] / 3e5 - 1 / 8) + abs(d101[j] / 3e5 - 1 / 8) +
                     abs(d011[j] / 3e5 - 1 / 8) + abs(d111[j] / 3e5 - 1 / 8))
  N.d[j] <- sum(ram.comb[[j]]$N.d) / 5e5
  N.u[j] <- sum(ram.comb[[j]]$N.u) / 5e5
  N.z[j] <- sum(ram.comb[[j]]$N.z) / 5e5
}

# Acceptance rate
mean(ram.comb.accept)

# N_discover
mean(rowSums(cbind(d100, d010, d001, d110, d101, d001) != 0))

# F_error
mean(avg.err)

# RAM's number of target density evaluations
colMeans(cbind(N.d, N.u, N.z))
sum(colMeans(cbind(N.d, N.u, N.z)))

######## Metropolis

metro.comb <- vector("list", 10)

for (j in 1 : 10) {
  metro.comb[[j]]  <- Metro(rep(0, test.dim), var.cov = var.c, 
                            n = 100, burn = 100, one.time.adaptation = TRUE)
}
# The sample size used in the article is 
# "sample.size = 1963200" and "burn.size = 1308800" for d = 3,
# "sample.size = 2261100" and "burn.size = 1507400" for d = 5,
# "sample.size = 2532300" and "burn.size = 1688200" for d = 7,
# "sample.size = 2840400" and "burn.size = 1893600" for d = 9, and
# "sample.size = 3210000" and "burn.size = 2140000" for d = 11.

# Summary
length.chain <- length(metro.comb[[1]]$accept)
metro.comb.accept <- rep(NA, 10)
d000 <- rep(NA, 10); d100 <- rep(NA, 10); d010 <- rep(NA, 10)
d001 <- rep(NA, 10); d110 <- rep(NA, 10); d101 <- rep(NA, 10)
d011 <- rep(NA, 10); d111 <- rep(NA, 10); avg.err <- rep(NA, 10)

for (j in 1 : 10) {
  metro.comb.accept[j] <- mean(metro.comb[[j]]$accept)
  d000[j] <- sum(metro.comb[[j]]$x[, 1] < 5 & metro.comb[[j]]$x[, 2] < 5 & metro.comb[[j]]$x[, 3] < 5)
  d100[j] <- sum(metro.comb[[j]]$x[, 1] > 5 & metro.comb[[j]]$x[, 2] < 5 & metro.comb[[j]]$x[, 3] < 5)
  d010[j] <- sum(metro.comb[[j]]$x[, 1] < 5 & metro.comb[[j]]$x[, 2] > 5 & metro.comb[[j]]$x[, 3] < 5)
  d001[j] <- sum(metro.comb[[j]]$x[, 1] < 5 & metro.comb[[j]]$x[, 2] < 5 & metro.comb[[j]]$x[, 3] > 5)
  d110[j] <- sum(metro.comb[[j]]$x[, 1] > 5 & metro.comb[[j]]$x[, 2] > 5 & metro.comb[[j]]$x[, 3] < 5)
  d101[j] <- sum(metro.comb[[j]]$x[, 1] > 5 & metro.comb[[j]]$x[, 2] < 5 & metro.comb[[j]]$x[, 3] > 5)
  d011[j] <- sum(metro.comb[[j]]$x[, 1] < 5 & metro.comb[[j]]$x[, 2] > 5 & metro.comb[[j]]$x[, 3] > 5)
  d111[j] <- sum(metro.comb[[j]]$x[, 1] > 5 & metro.comb[[j]]$x[, 2] > 5 & metro.comb[[j]]$x[, 3] > 5)
  avg.err[j] <- mean(abs(d000[j] / length.chain - 1 / 8) + abs(d100[j] / length.chain - 1 / 8) + 
                     abs(d010[j] / length.chain - 1 / 8) + abs(d001[j] / length.chain - 1 / 8) + 
                     abs(d110[j] / length.chain - 1 / 8) + abs(d101[j] / length.chain - 1 / 8) +
                     abs(d011[j] / length.chain - 1 / 8) + abs(d111[j] / length.chain - 1 / 8))
}

# Acceptance rate
mean(metro.comb.accept)

# N_discover
mean(rowSums(cbind(d100, d010, d001, d110, d101, d001) != 0))

# F_error
mean(avg.err)

######## Parallel tempering

temperature <- 2^(0:4)

pt.target <- function(temp, x) {
  log(target(x)) / temp
}

pt.comb <- vector("list", 10)

for (j in 1 : 10) {
  pt.comb[[j]] <- pt(initial.x = rep(0, test.dim), varcov = var.c, 
                     temp = temperature, swap.rate = 1, 
                     one.time.adaptation = TRUE,
                     n = 100, burn = 100)
}
# The sample size used in the article is 
# "sample.size = 280458" and "burn.size = 186971" for d = 3,
# "sample.size = 323014" and "burn.size = 215343" for d = 5,
# "sample.size = 361758" and "burn.size = 241171" for d = 7,
# "sample.size = 405772" and "burn.size = 270514" for d = 9, and
# "sample.size = 458572" and "burn.size = 305714" for d = 11,

# Summary
length.chain <- length(pt.comb[[1]]$accept[, 1])
pt.comb.accept <- rep(NA, 10)
d000 <- rep(NA, 10); d100 <- rep(NA, 10); d010 <- rep(NA, 10)
d001 <- rep(NA, 10); d110 <- rep(NA, 10); d101 <- rep(NA, 10)
d011 <- rep(NA, 10); d111 <- rep(NA, 10); avg.err <- rep(NA, 10)

for (j in 1 : 10) {
  pt.comb.accept[j] <- mean(pt.comb[[j]]$accept[, 1])
  d000[j] <- sum(pt.comb[[j]]$x[, 1] < 5 & pt.comb[[j]]$x[, 2] < 5 & pt.comb[[j]]$x[, 3] < 5)
  d100[j] <- sum(pt.comb[[j]]$x[, 1] > 5 & pt.comb[[j]]$x[, 2] < 5 & pt.comb[[j]]$x[, 3] < 5)
  d010[j] <- sum(pt.comb[[j]]$x[, 1] < 5 & pt.comb[[j]]$x[, 2] > 5 & pt.comb[[j]]$x[, 3] < 5)
  d001[j] <- sum(pt.comb[[j]]$x[, 1] < 5 & pt.comb[[j]]$x[, 2] < 5 & pt.comb[[j]]$x[, 3] > 5)
  d110[j] <- sum(pt.comb[[j]]$x[, 1] > 5 & pt.comb[[j]]$x[, 2] > 5 & pt.comb[[j]]$x[, 3] < 5)
  d101[j] <- sum(pt.comb[[j]]$x[, 1] > 5 & pt.comb[[j]]$x[, 2] < 5 & pt.comb[[j]]$x[, 3] > 5)
  d011[j] <- sum(pt.comb[[j]]$x[, 1] < 5 & pt.comb[[j]]$x[, 2] > 5 & pt.comb[[j]]$x[, 3] > 5)
  d111[j] <- sum(pt.comb[[j]]$x[, 1] > 5 & pt.comb[[j]]$x[, 2] > 5 & pt.comb[[j]]$x[, 3] > 5)
  avg.err[j] <- mean(abs(d000[j] / length.chain - 1 / 8) + abs(d100[j] / length.chain - 1 / 8) + 
                     abs(d010[j] / length.chain - 1 / 8) + abs(d001[j] / length.chain - 1 / 8) + 
                     abs(d110[j] / length.chain - 1 / 8) + abs(d101[j] / length.chain - 1 / 8) +
                     abs(d011[j] / length.chain - 1 / 8) + abs(d111[j] / length.chain - 1 / 8))
}

# Acceptance rate
mean(pt.comb.accept)

# N_discover
mean(rowSums(cbind(d100, d010, d001, d110, d101, d001) != 0))

# F_error
mean(avg.err)

############################ Example 3: Sensor location problem

######## Data

# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors. 
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)

# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0, 
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)

######## Target joint posterior density

norm2 <- function(loca, locb) {
  sqrt(sum((loca - locb)^2))
}

l.target <- function(loc, R = 0.3, sigma = 0.02, Ob, Os, Xb, Xs, Yb, Ys) {

  First.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Ob[j, i]) *
      (1 - exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Ob[j, i]) 
    })
    First.term <- c(First.term, TEMP)
  }

  Second.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                 loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Os[i, j]) *
      (1 - exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                      loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Os[i, j]) 
    })
    Second.term <- c(Second.term, TEMP)
  }

  First.obs.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      dnorm(Yb[j, i], mean = norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)]), 
                      sd = sigma)^Ob[j, i]
    })
    First.obs.term <- c(First.obs.term, TEMP)
  }

  Second.obs.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      dnorm(Ys[i, j], mean = norm2(loc[(2 * i -1) : (2 * i)], 
                                   loc[(2 * j -1) : (2 * j)]), 
                      sd = sigma)^Os[i, j]
    })
    Second.obs.term <- c(Second.obs.term, TEMP)
  }

  log.lik <- sum(log(c(First.term, Second.term, First.obs.term, Second.obs.term)))
  post <- log.lik + sum(dnorm(loc, mean = rep(0, 8), sd = rep(10, 8), log = TRUE))
  post

}

######## RAM

ram.kernel <- function(current.location, current.aux, loc.number, scale) {

  eps <- 10^(-308)
  accept <- 0
  x.c <- current.location 
  log.x.c.den <- l.target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  z.c <- current.aux
  log.z.c.den <- l.target(z.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.c.den <- exp(log.z.c.den)

  # downhill
  x.p1 <- x.c
  x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) {
    x.p1 <- x.c
    x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }

  # uphill
  x.p2 <- x.p1
  x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                   rnorm(2, 0, scale)
  log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) {
    x.p2 <- x.p1
    x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                     rnorm(2, 0, scale)
    log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }

  # downhill for N.d
  N.dz <- 1     # number of total downhill trials for estimate
  z <- x.p2
  z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                rnorm(2, 0, scale)
  log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.den <- exp(log.z.den)
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) {
    z <- x.p2
    z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, scale)
    log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    z.den <- exp(log.z.den)
    N.dz <- N.dz + 1
  }

  # accept or reject the proposal
  min.nu <- min(1, (x.c.den + eps) / (z.c.den + eps))
  min.de <- min(1, (x.p2.den + eps) / (z.den + eps))
  l.mh <- log.x.p2.den - log.x.c.den + log(min.nu) - log(min.de)

  if (l.mh > -rexp(1)) {
    x.c <- x.p2
    z.c <- z
    accept <- 1
  }

  c(x.c, z.c, N.d, N.u, N.dz, accept)
}

MHwG.RAM <- function(initial.loc, initial.aux, jump.scale, 
             Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 4)
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  Nz <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram.kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      aux.t <- TEMP[9 : 16]
      Nd[i, j] <- TEMP[17]
      Nu[i, j] <- TEMP[18]
      Nz[i, j] <- TEMP[19]
      accept[i, j] <- TEMP[20]
    }
    out[i, ] <- loc.t
  }
  print(Sys.time())
  list(x = out[-c(1 : n.burn), ], 
       accept = accept[-c(1 : n.burn), ],
       N.d = Nd[-c(1 : n.burn), ],
       N.u = Nu[-c(1 : n.burn), ],
       N.z = Nz[-c(1 : n.burn), ])

}

j.scale <- rep(1.08, 4)
system.time(res.ram <- MHwG.RAM(runif(8), runif(8), jump.scale = j.scale, 
                                Ob, Os, Xb, Xs, Yb, Ys, 
                                n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 200000" and "burn.size = 20000".

system.time(res.ram.den <- MHwG.RAM(runif(8), runif(8), jump.scale = j.scale, 
                                Ob, Os, Xb, Xs, Yb, Ys, 
                                n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 20000000" and "burn.size = 20000".

######## Metropolis

Metro.kernel <- function(current.location, loc.number, jump.scale) {

  accept <- 0
  x.p <- current.location
  x.p[(2 * loc.number - 1) : (2 * loc.number)] <- x.p[(2 * loc.number - 1) : (2 * loc.number)] + 
                                                  rnorm(2, 0, jump.scale)
  l.metro <- l.target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys) - 
             l.target(current.location, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)

  if (l.metro >  -rexp(1)) {
    current.location <- x.p
    accept <- 1
  }
  c(current.location, accept)

}

MHwG.Metro <- function(initial.loc, jump.scale, Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {

    for (j in 1 : 4) {
      TEMP <- Metro.kernel(loc.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      accept[i, j] <- TEMP[9]
    }

    out[i, ] <- loc.t

  }

  print(Sys.time())
  list(x = out[-c(1 : n.burn), ], accept = accept[-c(1 : n.burn), ])

}

j.scale <- rep(1.08, 4)
system.time(res.mt <- MHwG.Metro(initial.loc = runif(8), jump.scale = j.scale, 
                                 Ob, Os, Xb, Xs, Yb, Ys, n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 1967150" and "burn.size = 20000".

system.time(res.mt.den <- MHwG.Metro(initial.loc = runif(8), jump.scale = j.scale, 
                                 Ob, Os, Xb, Xs, Yb, Ys, n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 20000000" and "burn.size = 20000".

######## Tempered transitions

tt.kernel <- function(current.location, tempbase, nstep, loc.number, scale) {

    temperature <- tempbase^(1 : nstep)
    j.scale <- scale

    thresh.candi <- -rexp(nstep * 2)
    x.save <- matrix(NA, nrow = nstep * 2 + 1, ncol = 8)
    x.save[1, ] <- current.location
    tempupdown <- c(temperature[1 : nstep], temperature[nstep : 1])
    x.t <- current.location
    l.density.save <- rep(NA, nstep * 2 + 1)
    l.density.save[1] <- l.target(x.t, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)

    j.scale.serial <- c(j.scale * 1.2^{0 : (nstep - 1)}, j.scale * 1.2^{(nstep - 1) : 0})
    for (j in 1 : (2 * nstep)) {
      x.p <- x.t
      x.p[(2 * loc.number - 1) : (2 * loc.number)] <- 
         x.p[(2 * loc.number - 1) : (2 * loc.number)] + rnorm(2, mean = 0, sd = j.scale.serial[j])
      l.density.p <- l.target(x.p, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
      metro <- (l.density.p - l.density.save[j]) / tempupdown[j]
      if (metro > thresh.candi[j]) { 
        x.t <- x.p 
        l.density.save[j + 1] <- l.density.p
      } else {
        l.density.save[j + 1] <- l.density.save[j]
      }
      x.save[j + 1, ] <- x.t
    }

    E <- -l.density.save
    
    l.energy <- diff(1 / c(1, temperature)) %*% E[(2 * nstep + 1) : (nstep + 2)] -
                diff(1 / c(1, temperature)) %*% E[1 : nstep]  #LOG(exp(-(F_d-F_u)))

	if (l.energy > -rexp(1)) { 
	  current.location <- x.save[2 * nstep + 1, ] 
	}
	
	c(current.location)

}

MHwG.tt <- function(initial.loc, temp.base, n.step, jump.scale, 
                    Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {

  print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {
    for.accept <- loc.t
    for (j in 1 : 4) {
      TEMP <- tt.kernel(loc.t, temp.base, n.step, j, jump.scale)
      loc.t <- TEMP[1 : 8]
    }
    out[i, ] <- loc.t
    accept[i, ] <- 1 - (out[i, ] == for.accept)[seq(1, 7, by = 2)]
  }
  print(Sys.time())
  list(x = out[-c(1 : n.burn), ], accept = accept[-c(1 : n.burn), ])

}

j.scale <- 0.9
system.time(res.tt <- MHwG.tt(runif(8), 2, 3, jump.scale = j.scale, 
                              Ob, Os, Xb, Xs, Yb, Ys, n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 311192" and "burn.size = 20000".

system.time(res.tt.den <- MHwG.tt(runif(8), 2, 3, jump.scale = j.scale, 
                              Ob, Os, Xb, Xs, Yb, Ys, n.sample = 100, n.burn = 100))
# The sample size used in the article is 
# "sample.size = 311192" and "burn.size = 20000".

######## Summary

# Table 4

colMeans(res.ram$accept)
colMeans(res.ram$N.d)
colMeans(res.ram$N.u)
colMeans(res.ram$N.z)

colMeans(res.mt$accept)
colMeans(res.tt$accept)

# Figure 5

par(mfrow = c(4, 3), font = 2, font.lab = 2, font.axis = 2, cex = 1.2,
    mai = c(0.7, 0.9, 0.7, 0.3), mgp = c(2.5, 0.5, 0), las = 1)

plot(res.mt$x[, c(1, 2)], pch = 46, xlim = c(-0.2, 0.8), ylim = c(0.3, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Metropolis: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1.2)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)
plot(res.ram$x[, c(1, 2)], pch = 46, xlim = c(-0.2, 0.8), ylim = c(0.3, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1.2)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)
plot(res.tt$x[, c(1, 2)], pch = 46, xlim = c(-0.2, 0.8), ylim = c(0.3, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Tempered transitions: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1.2)
abline(v = 0.5748, lty = 2, lwd = 1)
abline(h = 0.9069, lty = 2, lwd = 1)

plot(res.mt$x[, c(3, 4)], pch = 46, xlim = c(-0.3, 1.3), ylim = c(-0.5, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Metropolis: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1.2)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)
plot(res.ram$x[, c(3, 4)], pch = 46, xlim = c(-0.3, 1.3), ylim = c(-0.5, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1.2)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)
plot(res.tt$x[, c(3, 4)], pch = 46, xlim = c(-0.3, 1.3), ylim = c(-0.5, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Tempered transitions: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1.2)
abline(v = 0.0991, lty = 2, lwd = 1)
abline(h = 0.3651, lty = 2, lwd = 1)

plot(res.mt$x[, c(5, 6)], pch = 46, xlim = c(0, 1), ylim = c(-0.1, 0.65),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Metropolis: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1.2)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)
plot(res.ram$x[, c(5, 6)], pch = 46, xlim = c(0, 1), ylim = c(-0.1, 0.65),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1.2)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)
plot(res.tt$x[, c(5, 6)], pch = 46, xlim = c(0, 1), ylim = c(-0.1, 0.65),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Tempered transitions: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1.2)
abline(v = 0.2578, lty = 2, lwd = 1)
abline(h = 0.1350, lty = 2, lwd = 1)

plot(res.mt$x[, c(7, 8)], pch = 46, xlim = c(-1.2, 1.9), ylim = c(-0.8, 1.9),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Metropolis: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1.2)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)
plot(res.ram$x[, c(7, 8)], pch = 46, xlim = c(-1, 2), ylim = c(-0.5, 2),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1.2)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)
plot(res.tt$x[, c(7, 8)], pch = 46, xlim = c(-1.2, 1.9), ylim = c(-0.8, 1.9),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("Tempered transitions: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1.2)
abline(v = 0.8546, lty = 2, lwd = 1)
abline(h = 0.0392, lty = 2, lwd = 1)

# Figure 6

par(mfrow = c(4, 3), font = 2, font.lab = 2, font.axis = 2, cex = 1.2,
    mai = c(0.7, 0.9, 0.7, 0.3), mgp = c(2.5, 0.5, 0), las = 1)

hist(res.mt$x[, 1], 50, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 11), xlim = c(-0.2, 0.75))
title(expression(bold(paste("Metropolis: ", x[11]))))
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
lines(density(res.mt.den$x[, 1]), lwd = 1.5)
abline(v = 0.5748, lty = 2, lwd = 2)

hist(res.ram$x[, 1], 50, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 11), xlim = c(-0.2, 0.75))
title(expression(bold(paste("RAM: ", x[11]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
index <- which(!is.na(res))
lines(density(res.ram.den$x[, 1]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.5748, lty = 2, lwd = 2)

hist(res.tt$x[, 1], 50, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 11), xlim = c(-0.2, 0.75))
title(expression(bold(paste("Tempered transition: ", x[11]))))
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
index <- which(!is.na(res))
lines(density(res.tt.den$x[, 1]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.5748, lty = 2, lwd = 2)

hist(res.mt$x[, 3], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 3), xlim = c(-0.15, 1.3))
index <- which(!is.na(res))
title(expression(bold(paste("Metropolis: ", x[21]))))
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.8, cex = 1.2, las = 0)
lines(density(res.mt.den$x[, 3]), lwd = 1.5)
abline(v = 0.1, lty = 2, lwd = 2)

hist(res.ram$x[, 3], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 3), xlim = c(-0.15, 1.3))
title(expression(bold(paste("RAM: ", x[21]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.8, cex = 1.2, las = 0)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
lines(density(res.ram.den$x[, 3]), lwd = 1.5)
abline(v = 0.1, lty = 2, lwd = 2)

hist(res.tt$x[, 3], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 3), xlim = c(-0.15, 1.3))
title(expression(bold(paste("Tempered transition: ", x[21]))))
mtext(side = 2, text = "Density", line = 1.8, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
index <- which(!is.na(res))
lines(density(res.tt.den$x[, 3]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.1, lty = 2, lwd = 2)

hist(res.mt$x[, 5], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 15), xlim = c(0.15, 0.85))
index <- which(!is.na(res))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
lines(density(res.mt.den$x[, 5]), lwd = 1.5)
title(expression(bold(paste("Metropolis: ", x[31]))))
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.2578, lty = 2, lwd = 2)

hist(res.ram$x[, 5], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 15), xlim = c(0.15, 0.85))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
index <- which(!is.na(res))
title(expression(bold(paste("RAM: ", x[31]))))
lines(density(res.ram.den$x[, 5]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.2578, lty = 2, lwd = 2)

hist(res.tt$x[, 5], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 15), xlim = c(0.15, 0.85))
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
index <- which(!is.na(res))
title(expression(bold(paste("Tempered transition: ", x[31]))))
lines(density(res.tt.den$x[, 5]), lwd = 1.5)
legend("top", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.2578, lty = 2, lwd = 2)

hist(res.mt$x[, 7], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 2.3), xlim = c(-1.1, 1.8))
index <- which(!is.na(res))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
lines(density(res.mt.den$x[, 7]), lwd = 1.5)
title(expression(bold(paste("Metropolis: ", x[41]))))
abline(v = 0.85, lty = 2, lwd = 2)

hist(res.ram$x[, 7], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 2.3), xlim = c(-1.1, 1.8))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[41]))))
lines(density(res.ram.den$x[, 7]), lwd = 1.5)
legend("topleft", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.85, lty = 2, lwd = 2)

hist(res.tt$x[, 7], 30, prob = TRUE, main = "", xlab = "", ylab = "",
     ylim = c(0, 2.3), xlim = c(-1.1, 1.8))
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
title(expression(bold(paste("Tempered transition: ", x[41]))))
lines(density(res.tt.den$x[, 7]), lwd = 1.5)
legend("topleft", c(expression(bold("Marginal dist.")), 
                    expression(bold("True location"))), lwd = 3, 
       bty = "n", seg.len = 0.9, cex = 0.9, lty = c(1, 3))
abline(v = 0.85, lty = 2, lwd = 2)

############################ Example 4: Strong lens time delay

post.theta <- function (data, X, delta, c, previous.theta, tau.jump, 
                        tau.thresh, tau.prior.a, tau.prior.b, sigma.prior.a, sigma.prior.b) {

  time <- data[, 1]
  time.d <- time - delta
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  time.diff <- diff(time.comb)
  leng.X <- length(X)

  mu <- previous.theta[1]
  sigma <- previous.theta[2]
  tau <- previous.theta[3]

  a.i <- exp( - time.diff / tau )   # i = 2, 3, ..., 2n
  leng.a <- length(a.i)  # 2n - 1

  # updating mu
  mu.mean <- ( X[1] + sum( (X[-1] - a.i * X[-leng.X]) / (1 + a.i) ) ) / 
                     ( 1 + sum( (1 - a.i) / (1 + a.i) ) )
  mu.sd <- sqrt( tau * sigma^2 / 2 /
                         ( 1 + sum( (1 - a.i) / (1 + a.i) ) ) )
  inv.cdf <- runif(1, min = pnorm(-30, mean = mu.mean, sd = mu.sd), 
                      max = pnorm(30, mean = mu.mean, sd = mu.sd))
  mu <- qnorm(inv.cdf, mean = mu.mean, sd = mu.sd)

  # updating sigma
  sigma <- sqrt( ( sigma.prior.b + (X[1] - mu)^2 / tau +  
                   sum( ( X[-1] - a.i * X[-leng.X] - mu * (1 - a.i) )^2 / 
                        (1 - a.i^2) ) / tau ) / 
                 rgamma(1, shape = length(time) + sigma.prior.a, scale = 1) )
 
  # updating tau
  log.post.tau <- function(t) {
    a.i.post <- exp( - time.diff / t )
    - (length(time) + 1 + tau.prior.a) * log(t) - 0.5 * sum( log(1 - a.i.post^2) ) - 
    ( tau.prior.b  + (X[1] - mu)^2 / sigma^2 +
      sum( ( X[-1] - a.i.post * X[-leng.X] - mu * (1 - a.i.post) )^2 / 
           (1 - a.i.post^2) ) / sigma^2 ) / t
  }

  tau.p <- exp(log(tau) + tau.jump)
  l.metrop <- log.post.tau(tau.p) - log.post.tau(tau)
  l.hastings <- log(tau.p) - log(tau)

  # Accept-reject
  if (l.metrop + l.hastings > tau.thresh) {
    tau <- tau.p
  }
   
  out <- c(mu, sigma, tau)  
  out
 
}

log.post.delta <- function(delta, data, theta, c, log, unif) {

  time <- data[, 1]
  leng.time <- length(time)

  if (delta < unif[1] | delta > unif[2]) {
    0
  } else if (theta[1] < -30 | theta[1] > 30) {
    0
  } else if (c <  -60 | c > 60) {
    0
  } else {

    lcA <- data[, 2]
    se.lcA <- data[, 3]
    lcB <- data[, 4]
    se.lcB <- data[, 5]

    if (log == TRUE) {
      # transform into magnitude scale 
      se.lcA <- se.lcA * 2.5 / lcA / log(10)
      lcA <- -2.5 * log(lcA, base = 10)
      se.lcB <- se.lcB * 2.5 / lcB / log(10)
      lcB <- -2.5 * log(lcB, base = 10)
      c.ini <- mean(lcB) - mean(lcA)
    }

    mu <- theta[1]
    sigma <- theta[2]
    tau <- theta[3]

    # sorting time given delta
    time.d <- time - delta
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    time.comb <- time.temp[ord]
    leng.time.comb <- length(time.comb)
  
    # indicator taking on 1 for X(t - delta) and 0 for X(t)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

    lc.temp <- c(lcA, lcB - c)
    lc.comb <- lc.temp[ord]
    se.lc.temp <- c(se.lcA, se.lcB)
    se.lc.comb <- se.lc.temp[ord]
   
    # x.star.i, i = 1, 2, ..., 2n
    x <- lc.comb - mu

    # omega.i, i = 1, 2, ..., 2n to be saved
    B <- rep(NA, leng.time.comb)

    # x.hat.i, i = 1, 2, ..., 2n to be saved
    mu.i <- rep(NA, leng.time.comb)
    mu.star.i <- rep(NA, leng.time.comb)

    # a.i, i = 2, ..., 2n
    a.i <- exp( -diff(time.comb) / tau)

    # omega.i, i = 1, 2, ..., 2n
    var0 <- tau * sigma^2 / 2
    B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0)

    for (k in 2 : leng.time.comb) {
      B[k] <- se.lc.comb[k]^2 / ( se.lc.comb[k]^2 + 
                                  a.i[k - 1]^2 * (1 - B[k - 1]) * se.lc.comb[k - 1]^2 +
                                  var0 * (1 - a.i[k - 1]^2) ) 
    }  

    # x.hat.i, i = 1, 2, ..., 2n
    mu.i[1] <- (1 - B[1]) * x[1]
    for (k in 2 : leng.time.comb) {
      mu.i[k] <- (1 - B[k]) * x[k] + B[k] * a.i[k - 1] * mu.i[k - 1]
    }

    mu.star.i[1] <- 0
    mu.star.i[2 : leng.time.comb] <- a.i * mu.i[-leng.time.comb]

    var.star.i <- se.lc.comb^2 / B

    # log-likelihood
    sum(dnorm(x, mean = mu.star.i, sd = sqrt(var.star.i), log = TRUE))

  }

}

post.X <- function(data, X, theta, delta, c, log) {

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]

  if (log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  mu <- theta[1]
  sigma <- theta[2]
  tau <- theta[3]

  # sorting time given delta
  time.d <- time - delta
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  time.comb <- time.temp[ord]
  leng.time.comb <- length(time.comb)
  
  # indicator taking on 1 for X(t - delta) and 0 for X(t)
  ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]

  lc.temp <- c(lcA, lcB - c)
  lc.comb <- lc.temp[ord]
  se.lc.temp <- c(se.lcA, se.lcB)
  se.lc.comb <- se.lc.temp[ord]
   
  # a.i, i = 1, ..., 2n
  time.comb.diff <- diff(time.comb)
  a.i <- exp( -time.comb.diff / tau)

  # x.i, i = 1, 2, ..., 2n
  X <- X - mu
  x <- lc.comb - mu

  # a.i, i = 2, ..., 2n
  time.comb.diff <- diff(time.comb)
  a.i <- exp( -time.comb.diff / tau )
  
  # B, i = 1, 2, ..., 2n to be saved
  B <- rep(NA, leng.time.comb)

  # mu.i, i = 1, 2, ..., 2n to be saved
  mu.i <- rep(NA, leng.time.comb)
  
  # shrinkages
  var0 <- tau * sigma^2 / 2
  B[1] <- se.lc.comb[1]^2 / (se.lc.comb[1]^2 + var0 * (1 - a.i[1]^2))
  mu.i[1] <- (1 - B[1]) * x[1] + B[1] * a.i[1] * X[2]
  X[1] <- rnorm(1, mean = mu.i[1], sd = sqrt(se.lc.comb[1]^2 * (1 - B[1])))

  for (k in 2 : (leng.time.comb - 1)) {
    B[k] <- se.lc.comb[k]^2 / 
              ( se.lc.comb[k]^2 + var0 * (1 - a.i[k - 1]^2) * (1 - a.i[k]^2) / (1 - (a.i[k - 1] * a.i[k])^2) )
    mu.i[k] <- (1 - B[k]) * x[k] + 
               B[k] * (a.i[k] * (1 - a.i[k - 1]^2) * X[k + 1] + a.i[k - 1] * (1 - a.i[k]^2) * X[k - 1]) / 
                      (1 - (a.i[k - 1] * a.i[k])^2)
    X[k] <- rnorm(1, mean = mu.i[k], sd = sqrt(se.lc.comb[k]^2 * (1 - B[k])))
  }

  B[leng.time.comb] <- se.lc.comb[leng.time.comb]^2 / 
                       (se.lc.comb[leng.time.comb]^2 + var0 * (1 - a.i[leng.time.comb - 1]^2))
  mu.i[leng.time.comb] <- (1 - B[leng.time.comb]) * x[leng.time.comb] + 
                          B[leng.time.comb] * a.i[leng.time.comb - 1] * X[leng.time.comb - 1]
  X[leng.time.comb] <- rnorm(1, mean = mu.i[leng.time.comb], 
                                sd = sqrt(se.lc.comb[leng.time.comb]^2 * (1 - B[leng.time.comb])))

  X + mu

}     

######## RAM

timedelay.ram <- function(data, theta.ini, delta.ini, delta.scale,
                          tau.jump.scale, tau.prior.shape, tau.prior.scale, 
                          sigma.prior.shape, sigma.prior.scale, Unif, 
                          Log = TRUE, 
                          adapt.delta = 0.05, adapt.tau = 0.05,
                          sample.size = 50, warmingup.size = 50) {

    ram.transition <- function(current.x, current.z, prop.scale, epsilon) {

    x.c <- current.x
    z.c <- current.z
    x.dim <- length(x.c)
    x.c.log.den <- log.post.delta(x.c, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)
    z.c.log.den <- log.post.delta(z.c, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)
    accept <- 0

    # downhill
    N.d <- 1    
    x.p1 <- rnorm(x.dim, x.c, prop.scale)
    x.p1.log.den <- log.post.delta(x.p1, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)

    while (runif(1) > (exp(x.c.log.den) + epsilon) / (exp(x.p1.log.den) + epsilon)) {
      x.p1 <- rnorm(x.dim, x.c, prop.scale)
      x.p1.log.den <- log.post.delta(x.p1, data, c(mu.t, sigma.t, tau.t),
                                 c.t, log = Log, unif = Unif)
      N.d <- N.d + 1
    }

    # uphill
    N.u <- 1
    x.p2 <- rnorm(x.dim, x.p1, prop.scale)
    x.p2.log.den <- log.post.delta(x.p2, data, c(mu.t, sigma.t, tau.t),
                                 c.t, log = Log, unif = Unif)

    while (runif(1) > (exp(x.p2.log.den) + epsilon) / (exp(x.p1.log.den) + epsilon)) {
      x.p2 <- rnorm(x.dim, x.p1, prop.scale)
      x.p2.log.den <- log.post.delta(x.p2, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)
      N.u <- N.u + 1
    }

    # downhill for N.d
    N.z <- 1
    z.p <- rnorm(x.dim, x.p2, prop.scale)
    z.p.log.den <- log.post.delta(z.p, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)
    while (runif(1) > (exp(x.p2.log.den) + epsilon) / (exp(z.p.log.den) + epsilon)) {
      z.p <- rnorm(x.dim, x.p2, prop.scale)
      z.p.log.den <- log.post.delta(z.p, data, c(mu.t, sigma.t, tau.t),
                               c.t, log = Log, unif = Unif)
      N.z <- N.z + 1
    }

    # accept or reject the proposal
    min.nu <- min(1, (exp(x.c.log.den) + epsilon) / (exp(z.c.log.den) + epsilon))
    min.de <- min(1, (exp(x.p2.log.den) + epsilon) / (exp(z.p.log.den) + epsilon))
    l.mh <- x.p2.log.den - x.c.log.den + log(min.nu) - log(min.de)

    if (l.mh > -rexp(1)) {
      x.c <- x.p2
      z.c <- z.p
      accept <- 1
    }

    c(x.c, z.c, accept, N.d, N.u, N.z)
  }

  print(Sys.time())
  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]
  c.ini <- mean(lcB) - mean(lcA)

  if (Log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini
  c.t <- c.ini  # c ini
  aux.t <- delta.ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  X.t <- c(lcA, lcB - c.t)[ord]  # X ini
  z.t <- delta.ini

  delta.out <- rep(NA, total.sample.size)
  N.d.out <- rep(NA, total.sample.size)
  N.u.out <- rep(NA, total.sample.size)
  N.z.out <- rep(NA, total.sample.size)

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  metro_mh <- rep(0, total.sample.size)
  
  tau.jumps <- tau.jump.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  

  delta.thresh <- -rexp(total.sample.size)
  delta.scale.adapt <- 1
  tau.scale.adapt <- 1
  epsilon <- 10^(-308)

  for (i in 1 : total.sample.size) {

    # delta and X(t) update
    temp <- ram.transition(delta.t, z.t, delta.scale, epsilon)

    delta.t <- delta.out[i] <- temp[1]
    z.t <- temp[2]
    delta.accept[i] <- temp[3]
    N.d.out[i] <- temp[4]
    N.u.out[i] <- temp[5]
    N.z.out[i] <- temp[6]

    if (delta.accept[i] == 1) { 
      X.t <- post.X(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
	                delta = delta.t, c = c.t, log = Log)
    }

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]

    # c update
    c.t.mean <- sum( (lc.comb - X.t) * ind / se.lc.comb^2 )  / sum( ind / se.lc.comb^2 )
    c.t.sd <- 1 / sqrt(sum( ind / se.lc.comb^2 ))
    inv.cdf <- runif(1, min = pnorm(-60, mean = c.t.mean, sd = c.t.sd), 
                        max = pnorm(60, mean = c.t.mean, sd = c.t.sd))
    c.t <- qnorm(inv.cdf, mean = c.t.mean, sd = c.t.sd)

    # theta update
    tau.jump.adapt <- tau.scale.adapt * tau.jumps[i]
    theta.update <- post.theta(data, X = X.t, delta = delta.t, c = c.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- theta.update[1]
    sigma.t <- theta.update[2]
    tau.t <- theta.update[3]

  }

  print(Sys.time())

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept[-c(1 : warmingup.size)]),
              delta.accept = mean(delta.accept[-c(1 : warmingup.size)]), 
              N.d = N.d.out[-c(1 : warmingup.size)], 
              N.u = N.u.out[-c(1 : warmingup.size)], 
              N.z = N.z.out[-c(1 : warmingup.size)])

  out

}

# set the place (working directory) where you put the data
setwd("/Users/hyungsuktak/Desktop/Dropbox/Data/Astrostat/")
dat <- read.csv("q0957usno.csv", header = TRUE)
rag <- max(dat[, 1]) - min(dat[, 1])
uniform <- c(-rag, rag)

delta.ini <- 0
delta.jump.scale <- 700

system.time(res.ram <- timedelay.ram(data = dat, theta.ini = c(mean(dat[, 2]), 0.01, 200),
                                    delta.ini = delta.ini, delta.scale = delta.jump.scale,
                                    tau.jump.scale = 1.75,
                                    tau.prior.shape = 1, tau.prior.scale = 1,
                                    sigma.prior.shape = 1, sigma.prior.scale = 2 * 10^-7,
                                    Unif = uniform, Log = FALSE,
                                    adapt.delta = 0, adapt.tau = 0,
                                    sample.size = 100, warmingup.size = 100))
# The sample size used in the article is 
# "sample.size = 7061612" and "burn.size = 50000".

######## Tempered transitions

timedelay.temper <- function(data, theta.ini, delta.ini, delta.scale,
                      tau.jump.scale, tau.prior.shape, tau.prior.scale, 
                      sigma.prior.shape, sigma.prior.scale, Unif, 
                      Log = TRUE, nstep = 5, tempbase = 8, 
                      sample.size = 1000, warmingup.size = 500) {

  print(Sys.time())
  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]
  c.ini <- mean(lcB) - mean(lcA)

  if (Log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini
  c.t <- c.ini  # beta.0 ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  X.t <- c(lcA, lcB - c.t)[ord]  # X ini

  delta.out <- rep(NA, total.sample.size)

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  
  tau.jumps <- tau.jump.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  

  delta.thresh <- -rexp(total.sample.size)
  tau.scale.adapt <- 1

  j.scale <- rep(delta.scale, nstep)
  for (i in 2 : nstep) {
    j.scale[i] <- j.scale[i - 1] * 1.2
  }
  delta.temp.scale <- c(j.scale, rev(j.scale))
  
  for (i in 1 : total.sample.size) {
    
    # delta and X(t) update
    temperature <- tempbase^(1 : nstep)
    candi <- rnorm(nstep * 2, sd = delta.temp.scale)
    thresh.candi <- -rexp(nstep * 2)
    delta.save <- rep(NA, nstep * 2 + 1)
    delta.save[1] <- delta.t
    density.save <- rep(NA, nstep * 2 + 1)
    density.save[1] <- log.post.delta(delta.t, data, c(mu.t, sigma.t, tau.t), 
                                      c.t, log = Log, unif = Unif)
    tempupdown <- c(temperature[1 : nstep], temperature[nstep : 1])
    for (j in 1 : (2 * nstep)) {
      delta.p <- delta.save[j] + candi[j]    
      density.p <- log.post.delta(delta.p, data, c(mu.t, sigma.t, tau.t), 
                                  c.t, log = Log, unif = Unif)
      metrop <- (density.p - density.save[j]) / tempupdown[j]
      if (metrop > thresh.candi[j]) { 
        delta.save[j + 1] <- delta.p 
        density.save[j + 1] <- density.p
      } else {
        delta.save[j + 1] <- delta.save[j]
        density.save[j + 1] <- density.save[j]
      }
    }

    E <- -density.save
    
    l.energy <- diff(1 / c(1, temperature)) %*% E[(2 * nstep + 1) : (nstep + 2)] - 
                diff(1 / c(1, temperature)) %*% E[1 : nstep]  #LOG(exp(-(F_d-F_u)))

    if (l.energy > delta.thresh[i]) { 
	  delta.t <- delta.save[2 * nstep + 1] 
	  delta.accept[i] <- 1
	  X.t <- post.X(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
	                delta = delta.t, c = c.t, log = Log)
    }
	
    delta.out[i] <- delta.t  

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]

    # c update
    c.t.mean <- sum( (lc.comb - X.t) * ind / se.lc.comb^2 )  / sum( ind / se.lc.comb^2 )
    c.t.sd <- 1 / sqrt(sum( ind / se.lc.comb^2 ))
    inv.cdf <- runif(1, min = pnorm(-56.74, mean = c.t.mean, sd = c.t.sd), 
                        max = pnorm(56.74, mean = c.t.mean, sd = c.t.sd))
    c.t <- qnorm(inv.cdf, mean = c.t.mean, sd = c.t.sd)
   
    tau.jump.adapt <- tau.scale.adapt * tau.jumps[i]

    # theta update
    theta.update <- post.theta(data, X = X.t, delta = delta.t, c = c.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- theta.update[1]
    sigma.t <- theta.update[2]
    tau.t <- theta.update[3]

  }

  print(Sys.time())

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept[-c(1 : warmingup.size)]),
              delta.accept.rate = mean(delta.accept[-c(1 : warmingup.size)]))

  out

}

# set the place (working directory) where you put the data
setwd("/Users/hyungsuktak/Desktop/Dropbox/Data/Astrostat/")
dat <- read.csv("q0957usno.csv", header = TRUE)
rag <- max(dat[, 1]) - min(dat[, 1])
uniform <- c(-rag, rag)
delta.ini <- 0
delta.jump.scale <- 500
nstep <- 5
tempbase <- 4

system.time(res.tt <- timedelay.temper(data = dat, theta.ini = c(mean(dat[, 2]), 0.01, 200),
                        delta.ini = delta.ini, delta.scale = delta.jump.scale,
                        tau.jump.scale = 1.75, tau.prior.shape = 1, tau.prior.scale = 1,
                        sigma.prior.shape = 1, sigma.prior.scale = 2 * 10^-7, Unif = uniform,
                        Log = FALSE, nstep = nstep, tempbase = tempbase,
                        sample.size = 100, warmingup.size = 100))
# The sample size used in the article is 
# "sample.size = 5000000" and "burn.size = 50000".

######## Metropolis

timedelay.mt <- function(data, theta.ini, delta.ini, delta.scale,
                         tau.jump.scale, tau.prior.shape, tau.prior.scale, 
                         sigma.prior.shape, sigma.prior.scale, Unif, 
                         Log = TRUE,
                         sample.size = 50, warmingup.size = 50) {

  print(Sys.time())
  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]
  c.ini <- mean(lcB) - mean(lcA)

  if (Log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini
  c.t <- c.ini  # c ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  X.t <- c(lcA, lcB - c.t)[ord]  # X ini

  delta.out <- rep(NA, total.sample.size)

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  
  tau.jumps <- tau.jump.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  

  delta.thresh <- -rexp(total.sample.size)
  delta.scale.adapt <- 1
  tau.scale.adapt <- 1

  for (i in 1 : total.sample.size) {

    # delta and X(t) update
	delta.p <- rnorm(1, mean = delta.t, sd = delta.scale)
	l.metrop <- log.post.delta(delta.p, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif) -
	            log.post.delta(delta.t, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif)
	if (l.metrop > delta.thresh[i]) { 
	  delta.t <- delta.p 
	  delta.accept[i] <- 1
	  X.t <- post.X(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
	                delta = delta.t, c = c.t, log = Log)
	}
	
	delta.out[i] <- delta.t  

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]

    # c update
    c.t.mean <- sum( (lc.comb - X.t) * ind / se.lc.comb^2 )  / sum( ind / se.lc.comb^2 )
    c.t.sd <- 1 / sqrt(sum( ind / se.lc.comb^2 ))
    inv.cdf <- runif(1, min = pnorm(-56.74, mean = c.t.mean, sd = c.t.sd), 
                        max = pnorm(56.74, mean = c.t.mean, sd = c.t.sd))
    c.t <- qnorm(inv.cdf, mean = c.t.mean, sd = c.t.sd)

    tau.jump.adapt <- tau.scale.adapt * tau.jumps[i]
    # theta update
    theta.update <- post.theta(data, X = X.t, delta = delta.t, c = c.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- theta.update[1]
    sigma.t <- theta.update[2]
    tau.t <- theta.update[3]

  }

  print(Sys.time())

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept[-c(1 : warmingup.size)]),
              delta.accept.rate = mean(delta.accept[-c(1 : warmingup.size)]))

  out

}

# set the place (working directory) where you put the data
setwd("/Users/hyungsuktak/Desktop/Dropbox/Data/Astrostat/")
rag <- max(dat[, 1]) - min(dat[, 1])
uniform <- c(-rag, rag)

delta.ini <- 0
dc <- 700

system.time(res.mt <- timedelay.mt(data = dat, theta.ini = c(mean(dat[, 2]), 0.01, 200),
                                    delta.ini = delta.ini, delta.scale = dc,
                                    tau.jump.scale = 1.7, 
                                    tau.prior.shape = 1, tau.prior.scale = 1, 
                                    sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7, 
                                    Unif = uniform, Log = FALSE, 
                                    sample.size = 100, warmingup.size = 100))
# The sample size used in the article is 
# "sample.size = 22216816" and "burn.size = 50000".

######## Metropolis with mixture jumping rule

timedelay.mt.mix <- function(data, theta.ini, delta.ini, delta.scale,
                      tau.jump.scale, tau.prior.shape, tau.prior.scale, 
                      sigma.prior.shape, sigma.prior.scale, Unif, 
                      Log = TRUE,
                      sample.size = 50, warmingup.size = 50) {

  print(Sys.time())
  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]
  c.ini <- mean(lcB) - mean(lcA)

  if (Log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini
  c.t <- c.ini  # c ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  X.t <- c(lcA, lcB - c.t)[ord]  # X ini

  delta.out <- rep(NA, total.sample.size)

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  
  tau.jumps <- tau.jump.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  

  delta.thresh <- -rexp(total.sample.size)
  delta.scale.adapt <- 1
  tau.scale.adapt <- 1

  for (i in 1 : total.sample.size) {

    # delta and X(t) update
    if (rbinom(1, 1, 0.5)) {
      delta.p <- rnorm(1, mean = delta.t, sd = delta.scale)
    } else {
      delta.p <- runif(1, Unif[1], Unif[2])
    }

	l.metrop <- log.post.delta(delta.p, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif) -
	            log.post.delta(delta.t, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif)
	if (l.metrop > delta.thresh[i]) { 
	  delta.t <- delta.p 
	  delta.accept[i] <- 1
	  X.t <- post.X(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
	                delta = delta.t, c = c.t, log = Log)
	}
	
	delta.out[i] <- delta.t  

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]

    # c update
    c.t.mean <- sum( (lc.comb - X.t) * ind / se.lc.comb^2 )  / sum( ind / se.lc.comb^2 )
    c.t.sd <- 1 / sqrt(sum( ind / se.lc.comb^2 ))
    inv.cdf <- runif(1, min = pnorm(-56.74, mean = c.t.mean, sd = c.t.sd), 
                        max = pnorm(56.74, mean = c.t.mean, sd = c.t.sd))
    c.t <- qnorm(inv.cdf, mean = c.t.mean, sd = c.t.sd)

    tau.jump.adapt <- tau.scale.adapt * tau.jumps[i]
    # theta update
    theta.update <- post.theta(data, X = X.t, delta = delta.t, c = c.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- theta.update[1]
    sigma.t <- theta.update[2]
    tau.t <- theta.update[3]

  }

  print(Sys.time())

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept[-c(1 : warmingup.size)]),
              delta.accept.rate = mean(delta.accept[-c(1 : warmingup.size)]))

  out

}

# set the place (working directory) where you put the data
setwd("/Users/hyungsuktak/Desktop/Dropbox/Data/Astrostat/")
dat <- read.csv("q0957usno.csv", header = TRUE)
rag <- max(dat[, 1]) - min(dat[, 1])
uniform <- c(-rag, rag)

delta.ini <- 0
delta.jump.scale <- 700

system.time(res.mt.mix <- timedelay.mt.mix(data = dat, theta.ini = c(mean(dat[, 2]), 0.01, 200),
                                    delta.ini = 0, delta.scale = delta.jump.scale,
                                    tau.jump.scale = 1.7, 
                                    tau.prior.shape = 1, tau.prior.scale = 1, 
                                    sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7, 
                                    Unif = uniform, Log = FALSE, 
                                    sample.size = 100, warmingup.size = 100))
# The sample size used in the article is 
# "sample.size = 19660188" and "burn.size = 50000".

######## An oracle M-H within Gibbs sampler 
######## (based on the information about the size and location of the mode.)

timedelay.oracle <- function(data, theta.ini, delta.ini,
                      tau.jump.scale, tau.prior.shape, tau.prior.scale, 
                      sigma.prior.shape, sigma.prior.scale, Unif, 
                      ASIS = FALSE, Log = TRUE,
                      sample.size = 50, warmingup.size = 50) {

  print(Sys.time())
  total.sample.size <- sample.size + warmingup.size

  time <- data[, 1]
  leng.time <- length(time)

  lcA <- data[, 2]
  se.lcA <- data[, 3]
  lcB <- data[, 4]
  se.lcB <- data[, 5]
  c.ini <- mean(lcB) - mean(lcA)

  if (Log == TRUE) {
    # transform into log scale
    se.lcA <- se.lcA * 2.5 / lcA / log(10)
    lcA <- -2.5 * log(lcA, base = 10)
    se.lcB <- se.lcB * 2.5 / lcB / log(10)
    lcB <- -2.5 * log(lcB, base = 10)
    c.ini <- mean(lcB) - mean(lcA)
  }

  delta.t <- delta.ini  # delta ini
  mu.t <- theta.ini[1]  # mu ini
  sigma.t <- theta.ini[2]  #sigma ini
  tau.t <- theta.ini[3]  #tau ini
  c.t <- c.ini  # c ini

  # sorting time given delta
  time.d <- time - delta.t
  time.temp <- c(time, time.d)
  ord <- order(time.temp)
  X.t <- c(lcA, lcB - c.t)[ord]  # X ini

  delta.out <- rep(NA, total.sample.size)

  tau.accept <- rep(0, total.sample.size)
  delta.accept <- rep(0, total.sample.size)
  
  tau.jumps <- tau.jump.scale * rnorm(total.sample.size)
  tau.thresh <- -rexp(total.sample.size)  

  delta.thresh <- -rexp(total.sample.size)
  tau.scale.adapt <- 1

  for (i in 1 : total.sample.size) {

    # delta and X(t) update
    if (rbinom(1, 1, 0.1)) {
      delta.p <- runif(1, 400, 450)
    } else {
      delta.p <- runif(1, 1050, Unif[2])
    }

	l.metrop <- log.post.delta(delta.p, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif) -
	            log.post.delta(delta.t, data, c(mu.t, sigma.t, tau.t), c.t, log = Log, unif = Unif)
	if (l.metrop > delta.thresh[i]) { 
	  delta.t <- delta.p 
	  delta.accept[i] <- 1
	  X.t <- post.X(data, X.t, theta = c(mu.t, sigma.t, tau.t), 
	                delta = delta.t, c = c.t, log = Log)
	}
	
	delta.out[i] <- delta.t  

    time.d <- time - delta.t
    time.temp <- c(time, time.d)
    ord <- order(time.temp)
    ind <- c(rep(0, leng.time), rep(1, leng.time))[ord]
    leng.X <- length(X.t)
    lc.comb <- c(lcA, lcB)[ord]
    se.lc.comb <- c(se.lcA, se.lcB)[ord]

    # c update
    c.t.mean <- sum( (lc.comb - X.t) * ind / se.lc.comb^2 )  / sum( ind / se.lc.comb^2 )
    c.t.sd <- 1 / sqrt(sum( ind / se.lc.comb^2 ))
    inv.cdf <- runif(1, min = pnorm(-56.74, mean = c.t.mean, sd = c.t.sd), 
                        max = pnorm(56.74, mean = c.t.mean, sd = c.t.sd))
    c.t <- qnorm(inv.cdf, mean = c.t.mean, sd = c.t.sd)

    if (ASIS == TRUE) {
      K.t <- X.t + c.t * ind
      K.t.cent <- K.t - mu.t
      time.comb <- time.temp[ord]
      time.diff <- diff(time.comb)
      # i = 2, 3, ..., 2n
      a.i <- exp( - time.diff / tau.t )

      y.c <- c(K.t.cent[1], K.t.cent[-1] - a.i * K.t.cent[-leng.X])
      X.c <- c(ind[1], ind[-1] - a.i * ind[-leng.X])
      V.inv.elem <- c(1, 1 / (1 - a.i^2))  #tau * sigma^2 /2 cancelled   
      XVX <- t(X.c * V.inv.elem) %*% X.c 
      c.mean <- t(X.c * V.inv.elem) %*% y.c  / XVX
      c.sd <- sqrt(tau.t * sigma.t^2 / 2 / XVX)
      inv.cdf <- runif(1, min = pnorm(-56.74, mean = c.mean, sd = c.sd), 
                          max = pnorm(56.74, mean = c.mean, sd = c.sd))
      c.t <- qnorm(inv.cdf, mean = c.mean, sd = c.sd)

      X.t <- K.t - c.t * ind    # synchronization
    }

    tau.jump.adapt <- tau.scale.adapt * tau.jumps[i]
    # theta update
    theta.update <- post.theta(data, X = X.t, delta = delta.t, c = c.t, 
                               previous.theta = c(mu.t, sigma.t, tau.t), 
                               tau.jump = tau.jump.adapt, tau.thresh = tau.thresh[i],
                               tau.prior.a = tau.prior.shape, tau.prior.b = tau.prior.scale, 
                               sigma.prior.a =  sigma.prior.shape, sigma.prior.b = sigma.prior.scale)

    if (theta.update[3] != tau.t) {
      tau.accept[i] <- 1
    }
 
    mu.t <- theta.update[1]
    sigma.t <- theta.update[2]
    tau.t <- theta.update[3]

  }

  print(Sys.time())

  out <- list(delta = delta.out[-c(1 : warmingup.size)],
              tau.accept.rate = mean(tau.accept[-c(1 : warmingup.size)]),
              delta.accept.rate = mean(delta.accept[-c(1 : warmingup.size)]))

  out

}

# set the place (working directory) where you put the data
setwd("/Users/hyungsuktak/Desktop/Dropbox/Data/Astrostat/")
dat <- read.csv("q0957usno.csv", header = TRUE)
rag <- max(dat[, 1]) - min(dat[, 1])
uniform <- c(-rag, rag)

delta.ini <- 0

system.time(res.oracle <- timedelay.oracle(data = dat, 
                                    theta.ini = c(mean(dat[, 2]), 0.01, 200),
                                    delta.ini = 0, 
                                    tau.jump.scale = 1.7, 
                                    tau.prior.shape = 1, tau.prior.scale = 1, 
                                    sigma.prior.shape = 1, sigma.prior.scale = 2 / 10^7, 
                                    Unif = uniform, Log = FALSE, 
                                    sample.size = 100, warmingup.size = 100))
# The sample size used in the article is 
# "sample.size = 20000000" and "burn.size = 50000".



######## Table 6

res.ram$delta.accept
res.tt$delta.accept.rate
res.mt$delta.accept.rate
res.mt.mix$delta.accept.rate

mean(res.ram$N.d)
mean(res.ram$N.u)
mean(res.ram$N.z)
mean(res.ram$N.z) + mean(res.ram$N.u) + mean(res.ram$N.d)

delta.middle <- res.ram$delta < 800
sum(diff(ifelse(delta.middle == "TRUE", 1, 0) == 1) != 0)

delta.middle <- res.tt$delta < 800
sum(diff(ifelse(delta.middle == "TRUE", 1, 0) == 1) != 0)

delta.middle <- res.mt$delta < 800
sum(diff(ifelse(delta.middle == "TRUE", 1, 0) == 1) != 0)

delta.middle <- res.mt.mix$delta < 800
sum(diff(ifelse(delta.middle == "TRUE", 1, 0) == 1) != 0)

######## Figure 8

par(mfrow = c(4, 1), font = 2, font.lab = 2, font.axis = 2, cex = 1.2,
    mai = c(0.7, 0.9, 0.7, 0.3), mgp = c(2.5, 0.5, 0), las = 1)

fig.height <- 0.05
hist(res.tt$delta, 3000, ylim = c(0, fig.height), prob = TRUE,
     main = "Tempered transitions", xlab = "")
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)

hist(res.ram$delta, 3000, ylim = c(0, fig.height), 
     main = "RAM", xlab = "", prob = TRUE)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)

hist(res.mt$delta, 3000, ylim = c(0, fig.height), 
     main = "Metropolis", xlab = "", prob = TRUE)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)

hist(res.mt.mix$delta, 3000, ylim = c(0, fig.height), 
     main = "Metropolis with mixture jumping rule", xlab = "", prob = TRUE)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)

# The right column

par(mfrow = c(4, 1), font = 2, font.lab = 2, font.axis = 2, cex = 1.2,
    mai = c(0.7, 0.9, 0.7, 0.3), mgp = c(2.5, 0.5, 0), las = 1)

hist(res.tt$delta, 3000, ylim = c(0, 0.05), xlim = c(420, 426),
     main = "Detail near 423 days", xlab = "", ylab = "", yaxt = "n", prob = TRUE)
lines(density(res.oracle$delta, n = 2^13, adjust = 0.1), lwd = 2, lty = 1)
mtext(side = 2, text = expression(bold("Density")), line = 2.3, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)
legend("topleft", c("Marginal posterior", "density"), lwd = 3, 
       col = c(1, 0), bty = "n", seg.len = 1, cex = 1, lty = c(1))
axis(2, at = seq(0, 0.05, by = 0.01), labels = TRUE)

hist(res.ram$delta, 3000, ylim = c(0, 0.05), xlim = c(420, 426),
     main = "", xlab = "", ylab = "", yaxt = "n", prob = TRUE)
lines(density(res.oracle$delta, n = 2^13, adjust = 0.1), lwd = 2, lty = 1)
mtext(side = 2, text = expression(bold("Density")), line = 2.3, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)
legend("topleft", c("Marginal posterior", "density"), lwd = 3, 
       col = c(1, 0), bty = "n", seg.len = 1, cex = 1, lty = c(1))
axis(2, at = seq(0, 0.05, by = 0.01), labels = TRUE)

hist(res.mt$delta, 3000, ylim = c(0, 0.05), xlim = c(420, 426),
     main = "", xlab = "", ylab = "", yaxt = "n", prob = TRUE)
lines(density(res.oracle$delta, n = 2^13, adjust = 0.1), lwd = 2, lty = 1)
mtext(side = 2, text = expression(bold("Density")), line = 2.3, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)
legend("topleft", c("Marginal posterior", "density"), lwd = 3, 
       col = c(1, 0), bty = "n", seg.len = 1, cex = 1, lty = c(1))
axis(2, at = seq(0, 0.05, by = 0.01), labels = TRUE)

hist(res.mt.mix$delta, 3000, ylim = c(0, 0.05), xlim = c(420, 426),
     main = "", xlab = "", ylab = "", yaxt = "n", prob = TRUE)
lines(density(res.oracle$delta, n = 2^13, adjust = 0.1), lwd = 2, lty = 1)
mtext(side = 2, text = expression(bold("Density")), line = 2.3, cex = 1.2, las = 0)
mtext(side = 1, text = expression(bold(Delta)), line = 1.5, cex = 1.2)
legend("topleft", c("Marginal posterior", "density"), lwd = 3, 
       col = c(1, 0), bty = "n", seg.len = 1, cex = 1, lty = c(1))
axis(2, at = seq(0, 0.05, by = 0.01), labels = TRUE)


