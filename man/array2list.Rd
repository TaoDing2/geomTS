% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{array2list}
\alias{array2list}
\title{Array object to list object}
\usage{
array2list(S)
}
\arguments{
\item{S}{An \eqn{p \times p \times n} array object.}
}
\value{
List object with length \eqn{n}. Note that each list contains a \eqn{p \times p} square matrix.
}
\description{
Transform array object to list object
}
\examples{
S = array(0,c(5,5,10))
for(i in 1:10){S[,,i] = CovM(5)}
X = array2list(S)

data(EEG)
S = EEG$seizure1
Slist = array2list(S)
length(Slist)
}
\seealso{
\code{\link[=list2array]{list2array()}}
}
