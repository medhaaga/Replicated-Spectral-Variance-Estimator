
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

