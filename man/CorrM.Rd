% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{CorrM}
\alias{CorrM}
\title{Full-rank correlation matrix}
\usage{
CorrM(p)
}
\arguments{
\item{p}{Row (or column) number of the symmetric matrix.}
}
\value{
A full-rank correlation matrix. It is also a SPD matrix.
}
\description{
A easy way to create an arbitrary \eqn{p \times p} correlation matrix with full rank from a covariance matrix.
}
\details{
It is computed by transforming a covariance matrix from \code{\link[=CovM]{CovM()}} into a correlation matrix.
}
\examples{
CorrM(5)
}
\seealso{
\code{\link[=CovM]{CovM()}}
}
