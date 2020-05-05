
## Implements the Markov chain Gibbs sampler
## for the Boomerang density
boom_gibbs <- function(nsim, start = c(0,0), A = 1, C = 7, B = 0)
{

    x1.samp <- rep(0, nsim)
    y1.samp <- rep(0, nsim)

    x1.samp[1] <- start[1]
    y1.samp <- start[2]

    for (i in 2:nsim)
    {
    x1.samp[i] <- rnorm(1,mean = C/(1 + A*y1.samp[i-1]^2), sd = sqrt(1/(1 + A*y1.samp[i-1]^2)))
    y1.samp[i] <- rnorm(1,mean = C/(1 + A*x1.samp[i]^2), sd = sqrt(1/(1 + A*x1.samp[i]^2)))
    }
    return(cbind(x1.samp, y1.samp))
}