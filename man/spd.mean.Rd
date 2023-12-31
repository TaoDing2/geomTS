% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AIfuncs.R
\name{spd.mean}
\alias{spd.mean}
\title{Frechet mean in affine invariant geometry}
\usage{
spd.mean(
  S,
  MaxIt = NULL,
  conver = NULL,
  dt = NULL,
  store.M = NULL,
  method = c("Sturm", "GD")
)
}
\arguments{
\item{S}{Array object with \eqn{p \times p \times n}.}

\item{MaxIt}{Maximum iterations of algorithm. The default value is 200.}

\item{conver}{Threshold of convergence. Default value is \eqn{10^{-4}}}

\item{dt}{Step size of algorithm. Normally, the default value is set  as 0.1.}

\item{store.M}{Logical values with \code{FALSE} as default value. If it is \code{TRUE}, we will store the estimated means and variance for each iteration, vice versa.}

\item{method}{Either "Sturm" (Sturm algorithm) nor "GD" (Gradient descent algorithm).
When it is \code{NULL}, gradient algorithm is chosen by the default.}
}
\value{
Estimated Frechet mean (and Frechet variance at each iteration if store.M = \code{TRUE})
}
\description{
Compute Frechet sample mean in affine invariant geometry \eqn{\mathcal{S}_p^+}
}
\examples{
S = lapply(1:10, function(i) CovM(5))
S = list2array(S)
FM = spd.mean(S,method = "Sturm",store.M = TRUE)
FM$mean

}
\references{
\insertRef{lenglet2006statistics}{geomTS}

\insertRef{sturm2003probability}{geomTS}
}
