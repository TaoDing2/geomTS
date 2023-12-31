% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QIfuncs.R
\name{cor.para_tran}
\alias{cor.para_tran}
\title{Parallel transport in quotient geometry}
\usage{
cor.para_tran(C1, C2, X1, D = NULL)
}
\arguments{
\item{C1}{A correlation matrix in \eqn{\mathcal{C}_p^+}.}

\item{C2}{A correlation matrix in \eqn{\mathcal{C}_p^+}.}

\item{X1}{A tangent vector matrix in \eqn{T_{C_1}\mathcal{C}_p^+}.
It is the initial velocity of horizontal geodesic at \eqn{\gamma^{\mathcal{H}}(t = 0) = C_1 \in \mathcal{S}_p^+}.}

\item{D}{A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
When it is \code{NULL}, we can obtain it from \code{\link[=D.optimal.position]{D.optimal.position()}}.}
}
\value{
A translated tangent vector matrix in \eqn{\mathcal{C}_p^+}.
}
\description{
Parallel transport the tangent vector \eqn{X \in T_{C_1}\mathcal{C}_p^+} into
another tangent space \eqn{T_{C_2}\mathcal{C}_p^+} along the geodesic connecting \eqn{C_1} and \eqn{C_2} in \eqn{\mathcal{C}_p^+}.
}
\details{
Given any two  points \eqn{C_1, C_2 \in \mathcal{C}_p^+} and tangent vector \eqn{X \in T_{C_1}\mathcal{C}_p^+},
we can  find the optimal position with respect to \eqn{C_1 \in \mathcal{S}_p^+} in the fibre of \eqn{C_2}
as \eqn{C_2^* = DC_2D \in \pi^{-1}(C_2)} which guarantees geodesic \eqn{\gamma(t), t \in [0,1]}
in \eqn{\mathcal{S}_p^+} joining \eqn{\gamma(0) = C_1 \in \mathcal{S}_p^+} and \eqn{\gamma(1) = C_2^* \in \mathcal{S}_p^+} is horizontal.
Henceforth, parallel transport of \eqn{X} from tangent space \eqn{T_{C_1}\mathcal{C}_p^+} to tangent space \eqn{T_{C_2}\mathcal{C}_p^+}
is defined as
\deqn{\mathcal{P}_{C_1 \rightarrow C_2}(X) = d\pi_{C_2^*} \left(  \mathcal{P}_{C_1 \rightarrow C_2^*}(\mathrm{Log}{C_1}(C_2^*)) \right)}
where \eqn{C_2^* \in \mathcal{S}_p^+} is the optimal position with respect to \eqn{C_1 \in \mathcal{S}_p^+}.
}
\note{
Parallel transport in quotient geometry in \eqn{\mathcal{C}_p^+} is uniquely determined.
}
\examples{
C1 = CorrM(5);C2 = CorrM(5)
X1 = cor.LogMap(C1,C2)
cor.para_tran(C1,C2,X1)
}
