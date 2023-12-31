% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE.R
\name{EstDiag}
\alias{EstDiag}
\title{Diagonal coefficients in the manifold-adapted model with diagonal coefficients}
\usage{
EstDiag(v, lagv)
}
\arguments{
\item{v}{A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.}

\item{lagv}{A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.}
}
\value{
Estimated diagonal coefficients, including \describe{
\item{Theta}{A list of estimated diagonal regressive coefficients.}
\item{Sigma}{A vector of estimated diagonals of standard deviation for the white noise}
}
}
\description{
Estimate diagonal coefficients in the manifold-adapted model with diagonal coefficients,
including coefficients in regressive terms and standard deviation of the white noise.
}
