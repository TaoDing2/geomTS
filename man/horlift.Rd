% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QIfuncs.R
\name{horlift}
\alias{horlift}
\title{Horizontal lift}
\usage{
horlift(X, S)
}
\arguments{
\item{X}{A tangent vector matrix in \eqn{T_C\mathcal{C}_p^+}.}

\item{S}{A SPD matrix in \eqn{\mathcal{S}_p^+}. It has the horizontal subspace \eqn{\mathcal{H}_S}.}
}
\value{
A horizontal symmetric matrix in horizontal subspace \eqn{\mathcal{H}_S}.
}
\description{
Horizontal lift the tangent vector \eqn{X \in T_C\mathcal{C}_p^+} to horizontal subspace \eqn{\mathcal{H}_S}.
}
\details{
Let \eqn{\pi: \mathcal{S}_p^+ \rightarrow \mathcal{C}_p^+} is a Riemannian submersion which is surjective onto \eqn{\mathcal{C}_p^+}.
Since \eqn{d\pi_S} is an isomorphism between \eqn{\mathcal{H}_S} and \eqn{T_C\mathcal{C}_p^+}, for every \eqn{C \in \mathcal{C}_p^+,
S \in \pi^{-1}(C) \in \mathcal{S}_p^+}, and \eqn{X \in T_C\mathcal{C}_p^+}, the unique horizontal lift \eqn{X^{\sharp}} in horizontal space
\eqn{\mathcal{H}_S} at \eqn{S \in \mathcal{S}_p^+} is defined as
\deqn{ X^{\sharp} = hor(\Delta_SX\Delta_S) \in \mathcal{H}_S .}
}
\examples{
S1 = CovM(5);S2 = CovM(5)
V1 = spd.LogMap(S1,S2)
X1 = deri_submer(V1,S1)
horlift(X1,S1)
}
