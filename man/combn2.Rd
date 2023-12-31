% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{combn2}
\alias{combn2}
\title{Squared matrix indices}
\usage{
combn2(p)
}
\arguments{
\item{p}{Integer}
}
\value{
A matrix containing indices of upper triangular squared matrix.
}
\description{
Generate squared matrix indices. It also generate all combinations of the elements of
\eqn{p} taken x at a time.
}
\examples{
combn2(5)
combn(5,2)
}
\seealso{
\code{\link[=combn]{combn()}}
}
