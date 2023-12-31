% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QIfuncs.R
\name{cor.geodesic}
\alias{cor.geodesic}
\title{Geodesic in quotient geometry}
\usage{
cor.geodesic(t, C1, C2, D = NULL)
}
\arguments{
\item{t}{Time. It should be in the interval of \eqn{[0,1]}.}

\item{C1}{A correlation matrix in \eqn{\mathcal{C}_p^+}. It is the starting point of geodesic \eqn{\c(t = 0) } in \eqn{\mathcal{C}_p^+}.}

\item{C2}{A correlation matrix in \eqn{\mathcal{C}_p^+}. It is the endpoint of geodesic \eqn{\c(t = 1) } in \eqn{\mathcal{C}_p^+}.}

\item{D}{A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
When it is \code{NULL}, we can obtain it from \code{\link[=D.optimal.position]{D.optimal.position()}}.}
}
\value{
A correlation matrix. It is in the geodesic \eqn{c(t)} at time \eqn{t, t \in [0,1]}
}
\description{
Compute geodesic \eqn{c(t), t \in [0,1]} in quotient geometry \eqn{\mathcal{C}_p^+}.
}
\details{
It is defined by \eqn{\mathrm{Exp}_{C_1}(t X)}, where \eqn{X = \mathrm{Log}_{C_1}(C_2)} and \eqn{t \in [0,1]}.
}
\note{
Any geodesic connecting two optimal positions is horizontal, and it is horizontal at one point will be horizontal at every point.
However, the points in segmentation of horizontal geodesic may not be in optimal positions but the tangent vector at base point are horizontal.
}
\examples{
C1 = CorrM(5); C2 = CorrM(5)
cor.geodesic(0.2,C1,C2)
}
