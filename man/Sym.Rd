% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{Sym}
\alias{Sym}
\title{Symmetrization operator}
\usage{
Sym(M)
}
\arguments{
\item{M}{A square matrix.}
}
\value{
A symmetric matrix.
}
\description{
Compute the Symmetrizatin operator:  \deqn{Sym(M) = \frac{1}{2}(M+M^T) }.
}
\examples{
M = matrix(c(1:9),3,3)
Sym(M)
}
