\name{combined_ccf}
\alias{combined_ccf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Replicated Cross Covariance/Correlation (R-CCF) Function Estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function `combined_ccf` estimates the CCF and R-CCF for a given list of parallel time series using the local mean and global mean respectively.
}
\usage{
combined_ccf <- function(x, center, chain, component, lag.max = NULL, type = c("covariance", "correlation"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a list of numeric matrices representing parallel time series
}
  \item{center}{
  global/overall arithmetic mean over all the replications
  }
  \item{chain}{
    numeric. The chain for which ACF and R-ACF are to be estimated
  }
  \item{component}{
  a numeric vector of two components for between which we want to calculate the CCF and R-CCF
  }
  \item{lag.max}{
  maximum lag at which to calculate the acf. Default is 10*log10(N) where N is the number of         observations. Will be automatically limited to one less than the number of observations in the serie.   }
  \item{type}{
  character string giving the type of acf to be computed. Allowed values are "correlation" (the default),and "covariance".
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function only calculates CCF and R-CCF between the given two components for a list of parallel time series.
For type = "correlation" and "covariance", the estimates are based on the sample covariance/correlation centered around sample mean of the concerned chain (for CCF) and overall sample mean (for R-CCF).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of two `acf` objects. First for CCF and  second for R-CCF.
}
\references{
%% ~put references to the literature/web site here ~
Our paper
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
