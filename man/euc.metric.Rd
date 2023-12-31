% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EUCfuncs.R
\name{euc.metric}
\alias{euc.metric}
\title{Euclidean distance in symmetric matrix space}
\usage{
euc.metric(S1, S2)
}
\arguments{
\item{S1}{A symmetric matrix.}

\item{S2}{A symmetric matrix.}
}
\value{
A non-negative number.
}
\description{
Calculate Euclidean distance for any two data points in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
}
\details{
Given any two points \eqn{S_1} and \eqn{S_2}  in symmetric matrix space \eqn{\mathrm{Sym}(p)}, we can define the Euclidean (or Frobenius) distance as
\deqn{d(S_1,S_2) =  \| S_1 -S_2 \|_F}
where \eqn{\|\cdot \|}  is the Frobenius distance.
}
\examples{
p = 5
S1 = CovM(p)
S2 = CovM(p)
euc.metric(S1,S2)
euc.metric(S2,S1)
}
