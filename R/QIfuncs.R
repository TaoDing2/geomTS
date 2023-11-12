################################################################################
###
###         Quotient geometry consisting of full-rank correlation matrix
###
################################################################################
#' Group action
#' @description Group action operation in SPD matrix space  \eqn{\mathcal{S}_p^+}.
#'
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}.
#' @param D A diagonal matrix with positive entries in Lie group \eqn{\mathcal{D}_p^+}.
#'
#' @return A SPD matrix after group action, which is also in manifold \eqn{\mathcal{S}_p^+}.
#' @export
#'
#' @details The Lie group \eqn{\mathcal{D}_p^+} of diagonal matrices with positive entries acts
#' smoothly, properly, and freely on Riemannian manifold \eqn{\mathcal{S}_p^+} via the action \eqn{\phi} at  point \eqn{S \in \mathcal{S}_p^+}:
#' \deqn{ \mathcal{D}_p^+ \times \mathcal{S}_p^+  \rightarrow \mathcal{S}_p^+}
#' i.e., \eqn{\phi_D(S) = DSD}, where \eqn{D \in \mathcal{D}_p^+}
#'
#' @references \insertRef{gallier2020differential}{geomTS}
#'
#' @examples S = CovM(5)
#' D = diag(c(1:5))
#' phi(S,D)
phi<- function(S,D){
  return(D%*%S%*%D)
}

#' Riemannian submersion
#' @description Riemannian submersion \eqn{\pi} between manifold \eqn{\mathcal{S}_p^+} and manifold \eqn{\mathcal{C}_p^+}.
#'
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}.
#'
#' @return A SPD matrix in quotient manifold \eqn{\mathcal{C}_p^+} which has 1 entries for all diagonals.
#' @export
#'
#' @details Let \eqn{\cdot :  \mathcal{D}_p^+ \times \mathcal{S}_p^+ \rightarrow \mathcal{S}_p^+} be a
#' smooth, free and proper action, with \eqn{ \mathcal{S}_p^+} a Lie group acting by isometries of \eqn{ \mathcal{S}_p^+}.
#' The smooth manifold \eqn{ \mathcal{C}_p^+ =  \mathcal{S}_p^+/ \mathcal{D}_p^+} can be seen as the quotient manifold via Riemannian submersion
#' \eqn{\pi:  \mathcal{S}_p^+ \rightarrow  \mathcal{C}_p^+}.
#' For example, a full rank correlation matrix from a covariance matrix is obtained by
#' \deqn{\pi : S \in \mathcal{S}_p^+ \rightarrow C =  \Delta_S^{-1}S\Delta_S^{-1} \in \mathcal{C}_p^+.}
#' where \eqn{\Delta_S = \sqrt{\mathrm{Diag}(S)}} and \eqn{\mathrm{Diag}(S) \in \mathcal{D}_p^+} is the diagonal matrix of \eqn{S}.
#'
#' @references \insertRef{david2019riemannian}{geomTS}
#'
#' @note In this package, we consider quotient manifold consisting of full-rank correlation matrices.
#' Therefore, correlation matrix [submer()] also could be achieved by function [cov2cor].
#'
#' @examples S = CovM(5)
#' submer(S)
submer <- function(S){
  D <- sol.sqrtm(Diag(S))
  return(phi(S,D))
}

#' Derivative of the Riemannian submersion
#' @description Differentiate Riemannian submersion \eqn{\pi} of \eqn{V \in T_S\mathcal{S}_p^+} at base point \eqn{S \in \mathcal{S}_p^+}.
#'
#' @param V A tangent vector matrix in \eqn{T_S\mathcal{S}_p^+}.
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}
#'
#' @return A tangent vector matrix in the tangent space \eqn{T_{\pi(S)}\mathcal{C}_p^+}.
#' @export
#'
#' @details The tangent vector matrix \eqn{X} is computed by the derivative of Riemannian submersion \eqn{d\pi_S(V)},
#' i.e. \eqn{X = d\pi_S(V)} is the  projected vector of \eqn{V} onto the tangent space \eqn{T_{\pi(S)}\mathcal{C}_p^+}.
#'
#' @note The tangent vector \eqn{X} by derivative of the Riemannian submersion \eqn{d\pi_S(V)} is a symmetric Hollow matrix with vanishing diagonals.
#'
#' @examples S1 = CovM(5)
#' S2 = CovM(5)
#' V1 = spd.LogMap(S1,S2)
#' deri_submer(V1,S1)
deri_submer <- function(V,S){
  D <- sol.sqrtm(Diag(S))
  tmp <- V - 0.5*(solve(Diag(S))%*%Diag(V)%*%S + S%*%Diag(V)%*%solve(Diag(S)))
  return(phi(tmp,D))
}

#' Vertical projection
#' @description The vertical component \eqn{W \in \mathcal{V}_S} of the tangent vector \eqn{V \in T_S\mathcal{S}_p^+}.
#'
#' @param V A tangent vector matrix in \eqn{T_S\mathcal{S}_p^+}.
#' @param S A SPD matrix in  \eqn{\mathcal{S}_p^+}.
#'
#' @return A symmetric matrix. It is the vertical component of tangent vector \eqn{V \in T_S\mathcal{S}_p^+}.
#' @export
#'
#' @details
#' Let \eqn{W \in \mathcal{V}_S} be the vertical component of the tangent vector \eqn{V \in T_S\mathcal{S}_p^+}. We have
#' \deqn{d\pi_S(W) = 0.}
#'
#' @note
#' The tangent sapce \eqn{T_{S}\mathcal{S_p^+}} to \eqn{\mathcal{S_p^+}} at \eqn{S} splits into the two components
#' \deqn{T_S\mathcal{S_p^+} =  \mathcal{H}_S \oplus  \mathcal{V}_S.}
#'
#' Let \eqn{W \in \mathcal{V}_S} be the vertical component of the tangent vector \eqn{V \in T_S\mathcal{S}_p^+}, and
#' \eqn{U \in \mathcal{H}_S} be the horizontal component of the tangent vector \eqn{V \in T_S\mathcal{S}_p^+}. We have
#' \deqn{V = U + W,}
#' and
#' \deqn{d\pi_S(W) = 0; d\pi_S(U) = d\pi_S(V) .}
#' This illustrates \eqn{d\pi_S} gives a linear isomorphism  between \eqn{\mathcal{H}_S} and \eqn{T_{\pi(S)}\mathcal{C}_p^+}.
#'
#' @seealso [hor()]
#' @examples S1 = CovM(5)
#' S2 = CovM(5)
#' V1 = spd.LogMap(S1,S2)
#' ver(V1,S1)
ver <- function(V,S){
  # V is the tangent vector in T_S at point S
  invS <- solve(S) # inverse of S
  I <- diag(dim(S)[1])     # identity matrix
  # calculate D
  # also compute as matrix formate as in thesis
  D = diag(c(solve(I + S*invS)%*%diag(invS%*%V)))
  verV <- S%*%D + D%*%S
  return(verV)
}

#' Horizontal projection
#' @description The horizontal component \eqn{U \in \mathcal{H}_S} of the tangent vector \eqn{V \in T_S\mathcal{S}_p^+}.
#'
#' @param V A tangent vector matrix in \eqn{T_S\mathcal{S}_p^+}.
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}.
#'
#' @return A symmetric matrix. It is the horizontal component of tangent vector \eqn{V \in T_S\mathcal{S}_p^+}.
#' @export
#'
#' @details
#' Let \eqn{U \in \mathcal{H}_S } be a vertical component of tangent vector \eqn{V \in T_S\mathcal{S}_p^+}. We have
#' \deqn{d\pi_S(V) = d\pi_S(U) = X,}
#' where \eqn{X} is the projected tangent vector in the tangent space \eqn{T_{\pi(S)}\mathcal{C}_p^+}.
#'
#' @note
#' Notes could be found in notes of [ver()].
#' @seealso [ver()]
#' @examples S1 = CovM(5); S2 = CovM(5)
#' V1 = spd.LogMap(S1,S2)
#' hor(V1,S1)
hor <- function(V,S){
  # Vertical vector
  verV <- ver(V,S)
  return(V - verV)
}

#'  Horizontal lift
#' @description
#' Horizontal lift the tangent vector \eqn{X \in T_C\mathcal{C}_p^+} to horizontal subspace \eqn{\mathcal{H}_S}.
#'
#' @param X A tangent vector matrix in \eqn{T_C\mathcal{C}_p^+}.
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}. It has the horizontal subspace \eqn{\mathcal{H}_S}.
#'
#' @return A horizontal symmetric matrix in horizontal subspace \eqn{\mathcal{H}_S}.
#' @export
#'
#' @details
#' Let \eqn{\pi: \mathcal{S}_p^+ \rightarrow \mathcal{C}_p^+} is a Riemannian submersion which is surjective onto \eqn{\mathcal{C}_p^+}.
#' Since \eqn{d\pi_S} is an isomorphism between \eqn{\mathcal{H}_S} and \eqn{T_C\mathcal{C}_p^+}, for every \eqn{C \in \mathcal{C}_p^+,
#' S \in \pi^{-1}(C) \in \mathcal{S}_p^+}, and \eqn{X \in T_C\mathcal{C}_p^+}, the unique horizontal lift \eqn{X^{\sharp}} in horizontal space
#' \eqn{\mathcal{H}_S} at \eqn{S \in \mathcal{S}_p^+} is defined as
#' \deqn{ X^{\sharp} = hor(\Delta_SX\Delta_S) \in \mathcal{H}_S .}
#'
#' @examples S1 = CovM(5);S2 = CovM(5)
#' V1 = spd.LogMap(S1,S2)
#' X1 = deri_submer(V1,S1)
#' horlift(X1,S1)
horlift <- function(X,S){
  D <- sqrtm(Diag(S))
  # horizontal vectors
  return(hor(phi(X,D),S))
}

#' Find the optimal position
#' @description
#' Find the optimal position for any two points \eqn{S_1 \in \pi^{-1}(C_1), S_2 \in \pi^{-1}(C_2)}
#' in \eqn{\mathcal{S}_p^+} via greedy descent algorithm.
#'
#' @param C1 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param C2 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param MaxIt Maximum iterations. Default is 200.
#' @param delta Step size. Default is 0.001.
#' @param conver Convergence parameter. Default is 4, which means the distance between two iterations is less then \eqn{10^{-4}}.
#' @param store.D Logical values. Default is \code{FALSE} which means we don't want to store the iterated \eqn{D} and the values of objective function.
#'
#' @return A diagonal matrix (iterated \eqn{D} and values of objective function if store.M = \code{TRUE}).
#' @export
#'
#' @details Suppose \eqn{C_1} and \eqn{C_2} both in \eqn{\mathcal{C}_p^+},
#' \eqn{C_1} lies in \eqn{\mathcal{S}_p^+} as it is a SPD matrix.
#' We aim to find a matrix \eqn{D \in \mathcal{D}_p^+}
#' which makes affine invariant metric of \eqn{C_1} and \eqn{DC_2D} in \eqn{\mathcal{S}_p^+} shortest.
D.optimal.position <- function(C1,C2,MaxIt= NULL,delta = NULL,conver = NULL,store.D = NULL){
  if(is.null(delta)) delta <- 0.1
  if(is.null(store.D)) store.D <- F
  if(is.null(conver)) conver <- 4
  if(is.null(MaxIt)) MaxIt <- 200
  # Initial D0 as identity matrix
  p = ncol(C1)
  D0 <- diag(p)
  # store diagonal matrices
  if(store.D){
    itD <- list()
  }
  # store distance
  d <- c()
  # iteration
  k <- 1
  convergence <- 1
  while (convergence > 10^{-conver}) {
    # Delta_t
    Delta  <- Diag(2*Sym(D0%*%logm(C2%*%D0%*%solve(C1)%*%D0)))
    # update diagonal matrix D
    D <- phi(expm(-delta*phi(Delta,sol.sqrtm(D0))),sqrtm(D0))
    if(store.D){
      itD[[k]] <- D
    }
    ### check the converegnce
    d[k] <- spd.metric(C1,phi(C2,D))
    if(k > 1) {
      convergence <- abs(d[k]- d[k-1])
    }
    ### Updating
    k = k + 1
    D0 = D
    ### break the loop when iteration reach the maximum iterations
    if(k > MaxIt) break
  }
  # find D w.r.t the optimal position
  if(store.D) {
    results = list(D = D0, itD = itD, dist = d)
  } else {
    results = D0
  }
  return(results)
}

#' Geodesic distance in  quotient geometry
#' @description
#' Given any two points \eqn{C_1} and \eqn{C_2} in quotient geometry \eqn{\mathcal{C}_p^+},
#' we compute the geodesic distance between these two points based on affine invariant metric in \eqn{\mathcal{S}_p^+}.
#'
#' @param C1 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param C2 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()]
#'
#' @details
#' Given any two points \eqn{C_1,C_2 \in \mathcal{C}_p^+} and chosen \eqn{\pi^{-1}(C_1) = I_pC_1I_p}, the quotient metric in quotient manifold is defined a
#' \deqn{d^{\mathrm{quo}}(C_1,C_2) = \inf \limits_{D \in \mathcal{D}_p^+} d^{\mathrm{ai}}(C_1,DC_2D) }
#'
#' @return A numerical value of distance between \eqn{C_1} and \eqn{C_2} in \eqn{\mathcal{C}_p^+}
#' @export
#'
#' @examples C1 = CorrM(5) ; C2 = CorrM(5)
#' cor.metric(C1,C2)
cor.metric <- function(C1,C2,D = NULL){
  # optimal position to C1
  if(is.null(D)) D <- D.optimal.position(C1,C2)
  S2 = phi(C2,D)
  return(spd.metric(C1,S2))
}

#' Inner product in quotient geometry
#' @description
#' Compute the inner product for any two tangent vector matrices \eqn{X, Y \in T_C\mathcal{C}_p^+} at the base point \eqn{C \in \mathcal{C}_p^+ }.
#'
#' @param X A tangent vector matrix in \eqn{T_{C}\mathcal{C}_p^+}.
#' @param Y A tangent vector matrix in \eqn{T_{C}\mathcal{C}_p^+}.
#' @param S A SPD matrix. It lies in the manifold \eqn{\mathcal{S}_p^+} and has horizontal subspace \eqn{\mathcal{H}_S}.
#' Also, \eqn{C = \pi(S) \in \mathcal{C}_p^+} and \eqn{C \in \mathcal{C}_p^+.}
#'
#' @return A numerical distance.
#' @export
#' @details For any \eqn{C \in \mathcal{C}_p^+}, \eqn{X,Y \in T_C\mathcal{C}_p^+}, and \eqn{S \in \pi^{-1}(C) \subset \mathcal{S}_p^+},
#' there exists unique horizontal lift \eqn{X^{\sharp}} and \eqn{Y^{\sharp}}  both in horizontal subspace
#' \eqn{\mathcal{H}_S} such that \eqn{d\pi_S(X^{\sharp}) = X} and \eqn{d\pi_S(Y^{\sharp}) = Y}.
#'
#' The quotient metric in quotient manifold \eqn{\mathcal{C}_p^+} is defined based on the affine invariant metric in \eqn{\mathcal{S}_p^+ } as
#' \deqn{ g^{\mathrm{quo}}_C\langle  X,Y\rangle  = g^{\mathrm{ai}}_S\langle  X^{\sharp},Y^{\sharp} \rangle .}
#' @note \itemize{
#' \item Two tangent vectors \eqn{X,Y} should both lie in the same tangent space \eqn{T_{C}\mathcal{C}_p^+}
#' \item Quotient metric \eqn{g^{\mathrm{quo}}_C} does not depend on the choice of \eqn{S} in the fiber of \eqn{\pi^{-1}(C) \subset \mathcal{S}_p^+}.
#' }
#'
#' @examples
#' C1 = CorrM(5); C2 = CorrM(5); C3 = CorrM(5)
#' X = cor.LogMap(C1,C2) ;Y = cor.LogMap(C1,C3)
#' cor.inner_product(X,Y,C1)
cor.inner_product <- function(X,Y,S){
  D <- sqrtm(Diag(S))
  X_S <- phi(X,D)
  Y_S <- phi(Y,D)
  verX <- ver(X_S,S)
  verY <- ver(Y_S,S)
  d <- spd.inner_product(X_S,Y_S,S) - spd.inner_product(verX,verY,S)
  return(d)
}


#' Logarithm maps in quotient geometry
#' @description Logarithm map in quotient geometry \eqn{\mathcal{C}_p^+}.
#' It aims to obtain a tangent vector matrix by computing the velocity from point \eqn{C_1} to point \eqn{C_2} both locating in \eqn{\mathcal{C}_p^+}.
#'
#' @param C1 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param C2 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()].
#'
#' @return A tangent vector matrix in \eqn{T_{C_1}\mathcal{C}_p^+}.
#'
#' @details Given any \eqn{C_1,C_2 \in \mathcal{C}_p^+}, we can  find the optimal position  in the fibre of \eqn{\pi^{-1}(C_2) \subset \mathcal{S}_p^+}
#' with respect to \eqn{C_1 \in \mathcal{S}_p^+} by optimizing objective function,
#' i.e., obtaining \eqn{C_2^* = DC_2D \in \pi^{-1}(C_2)  \subset \mathcal{S}_p^+},
#' so that \eqn{\mathrm{Log}_{C_1}(C_2^*) \in T_{C_1}\mathcal{C}_p^+} is a tangent vector,
#' also a horizontal vector, i.e. \eqn{ \mathrm{Log}{C_1}(C_2^*) = hor(\mathrm{Log}_{C_1}(C_2^*))}.
#' Consequently, we can define the logarithmic map in quotient manifold as
#' \deqn{\mathrm{Log}_{C_1}(C_2) =  d\pi_{C_1}(hor(\mathrm{Log}_{C_1}(C^*_2))) = d\pi_{C_1}(\mathrm{Log}_{C_1}(C^*_2)) \in T_{C_1}\mathcal{C}_p^+ }
#'
#' @note \itemize{
#' \item The matrix in \eqn{T_{C}\mathcal{C}_p^+} is a symmetric Hollow matrix with vanishing diagonals.
#' \item The logarithm map  \eqn{X = \mathrm{Log}_{C_1}(C_2)} in quotient manifold is uniquely determined.
#' \item The tangent vector is the initial velocity at the starting point in the geodesic connecting \eqn{C_1} and \eqn{C_2} both in \eqn{\mathcal{C}_p^+}.
#' \item It is the inverse map of [cor.ExpMap()].
#' }
#'
#' @seealso [cor.ExpMap()]
#' @export
#' @examples
#'  C1 = CorrM(5) ; C2 = CorrM(5)
#'  cor.LogMap(C1,C2)
#'  D = D.optimal.position(C1,C2)
#'  cor.LogMap(C1,C2,D)
cor.LogMap <- function(C1,C2,D = NULL){
  # optimal position w.r.t C1
  if(is.null(D)) D <- D.optimal.position(C1,C2)
  # tangent vector in horizontal space at C1
  V <- spd.LogMap(C1,phi(C2,D))
  # derivative of submersion
  return(deri_submer(V,C1))
}

#' Exponential map in quotient geometry
#' @description Exponential map in quotient geometry \eqn{\mathcal{C}_p^+}.
#' It projects the tangent vector matrix \eqn{X \in T_{C}\mathcal{C}_p^+} at point \eqn{C \in \mathcal{C}_p^+} back to \eqn{\mathcal{C}_p^+}.
#'
#' @param C A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param X A tangent vector matrix in \eqn{T_{C}\mathcal{C}_p^+}.
#'
#' @return A correlation matrix. It is the mapped point of \eqn{X} at base point \eqn{C} in the manifold \eqn{\mathcal{C}_p^+}.
#' @details For any point \eqn{C} in quotient manifold \eqn{\mathcal{C}_p^+} and tangent vector \eqn{X \in T_C\mathcal{C}_p^+} at the base point \eqn{C},
#'  we can "lift" \eqn{X} from \eqn{T_C\mathcal{C}_p^+} to horizontal space \eqn{\mathcal{H}_C} and then do exponential map in \eqn{\mathcal{S}_p^+} .
#'  Finally, Riemannian submersion is used to project this point onto quotient manifold so that exponential map in \eqn{T_C\mathcal{C}_p^+} is defined as
#'  \deqn{\mathrm{Exp}_{C}(X) = \pi(\mathrm{Exp}_{C}(X^{\sharp})) = \pi (C^{1/2}\exp (C^{-1/2}hor(X) C^{-1/2})C^{1/2}) }
#'
#' @note It is the inverse map of [cor.LogMap()]
#' @export
#' @examples
#' C1 = CorrM(5); C2 = CorrM(5)
#' X1 = cor.LogMap(C1,C2)
#' cor.ExpMap(C1,X1)
cor.ExpMap <- function(C,X){
  #### horizontal lift to horizontal space H_C
  horV <- horlift(X,C)
  #### exponential map in SPD
  S <- spd.ExpMap(C,horV)
  ### submersion onto COR
  return(submer(S))
}

#' Geodesic in quotient geometry
#' @description
#' Compute geodesic \eqn{c(t), t \in [0,1]} in quotient geometry \eqn{\mathcal{C}_p^+}.
#'
#' @param t Time. It should be in the interval of \eqn{[0,1]}.
#' @param C1 A correlation matrix in \eqn{\mathcal{C}_p^+}. It is the starting point of geodesic \eqn{\c(t = 0) } in \eqn{\mathcal{C}_p^+}.
#' @param C2 A correlation matrix in \eqn{\mathcal{C}_p^+}. It is the endpoint of geodesic \eqn{\c(t = 1) } in \eqn{\mathcal{C}_p^+}.
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()].
#'
#' @return A correlation matrix. It is in the geodesic \eqn{c(t)} at time \eqn{t, t \in [0,1]}
#' @export
#' @details It is defined by \eqn{\mathrm{Exp}_{C_1}(t X)}, where \eqn{X = \mathrm{Log}_{C_1}(C_2)} and \eqn{t \in [0,1]}.
#'
#' @note Any geodesic connecting two optimal positions is horizontal, and it is horizontal at one point will be horizontal at every point.
#' However, the points in segmentation of horizontal geodesic may not be in optimal positions but the tangent vector at base point are horizontal.
#'
#' @examples C1 = CorrM(5); C2 = CorrM(5)
#' cor.geodesic(0.2,C1,C2)
cor.geodesic <- function(t,C1,C2,D =NULL){
  # tangent vector at T_C0
  if(is.null(D)) D <- D.optimal.position(C1,C2)
  S2 = phi(C2,D)
  # geodesic in S points in optimal positions
  geo <- spd.geodesic(t,C1,S2)
  ### return the submersion
  return(submer(geo))
}

#' Parallel transport in quotient geometry
#' @description
#' Parallel transport the tangent vector \eqn{X \in T_{C_1}\mathcal{C}_p^+} into
#' another tangent space \eqn{T_{C_2}\mathcal{C}_p^+} along the geodesic connecting \eqn{C_1} and \eqn{C_2} in \eqn{\mathcal{C}_p^+}.
#'
#' @param C1 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param C2 A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param X1 A tangent vector matrix in \eqn{T_{C_1}\mathcal{C}_p^+}.
#' It is the initial velocity of horizontal geodesic at \eqn{\gamma^{\mathcal{H}}(t = 0) = C_1 \in \mathcal{S}_p^+}.
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()].
#'
#' @return A translated tangent vector matrix in \eqn{\mathcal{C}_p^+}.
#' @details
#' Given any two  points \eqn{C_1, C_2 \in \mathcal{C}_p^+} and tangent vector \eqn{X \in T_{C_1}\mathcal{C}_p^+},
#' we can  find the optimal position with respect to \eqn{C_1 \in \mathcal{S}_p^+} in the fibre of \eqn{C_2}
#'  as \eqn{C_2^* = DC_2D \in \pi^{-1}(C_2)} which guarantees geodesic \eqn{\gamma(t), t \in [0,1]}
#'   in \eqn{\mathcal{S}_p^+} joining \eqn{\gamma(0) = C_1 \in \mathcal{S}_p^+} and \eqn{\gamma(1) = C_2^* \in \mathcal{S}_p^+} is horizontal.
#'   Henceforth, parallel transport of \eqn{X} from tangent space \eqn{T_{C_1}\mathcal{C}_p^+} to tangent space \eqn{T_{C_2}\mathcal{C}_p^+}
#'   is defined as
#'   \deqn{\mathcal{P}_{C_1 \rightarrow C_2}(X) = d\pi_{C_2^*} \left(  \mathcal{P}_{C_1 \rightarrow C_2^*}(\mathrm{Log}{C_1}(C_2^*)) \right)}
#'   where \eqn{C_2^* \in \mathcal{S}_p^+} is the optimal position with respect to \eqn{C_1 \in \mathcal{S}_p^+}.
#' @note Parallel transport in quotient geometry in \eqn{\mathcal{C}_p^+} is uniquely determined.
#' @export
#'
#' @examples
#' C1 = CorrM(5);C2 = CorrM(5)
#' X1 = cor.LogMap(C1,C2)
#' cor.para_tran(C1,C2,X1)
cor.para_tran <- function(C1,C2,X1,D = NULL){
  # optimal position to C0
  if(is.null(D)) D <- D.optimal.position(C1,C2)
  ### horizontal lift
  horV <- horlift(X1,C1)
  ### optimal position of C2 w.r.t C0 in S
  S2 <- phi(C2,D)
  ### parallel transport in horizontal space
  ParTra <- spd.para_tran(C1,S2,horV)
  return(deri_submer(ParTra,S2))
}

#' Orthonormal basis at any base point in quotient geometry
#' @description
#' Construct orthonormal basis in the tangent space \eqn{T_{C}\mathcal{C}_p^+}.
#'
#' @param C A correlation matrix in \eqn{\mathcal{C}_p^+}
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()].
#'
#' @return A list of orthonormal basis in the tangent space \eqn{T_{C}\mathcal{C}_p^+}
#' containing \eqn{m = p(p - 1)/2} observations of \eqn{p \times p} basis.
#' @details
#' Let  \eqn{C} be any point in quotient manifold \eqn{\mathcal{C}_p^+}. The orthonormal basis in tangent space at base point \eqn{C} is defined as
#' \deqn{ \begin{aligned}
#'  \boldsymbol{E}^C  & = \mathcal{P}_{I_p \rightarrow C}(\boldsymbol{E}^{I_p}) \\
#'  & =  d \pi_{C^{\dag}}\left[(C^{\dag} I_p^{-1})^{1/2}(\boldsymbol{E}^{I_p})^{\sharp} ((C^{\dag} I_p^{-1})^{1/2})^T \right]\\
#'  & = d \pi_{C^{\dag}}\left[{C^{\dag}}^{1/2}hor(\boldsymbol{E}) {C}^{1/2} \right]
#'  \end{aligned}
#'  }
#' where \eqn{\boldsymbol{E}^{I_p}} is the orthonormal basis in the tangent space \eqn{T_{I_p}\mathcal{C}_p^+} with the dimension of \eqn{m = p(p-1)/2},
#' and it is equivalent to orthonormal basis in horizontal space \eqn{\mathcal{H}_{I_p}} in [hor.ortho_basis()].
#' And, \eqn{C^{\dag} = DCD} is the optimal position with respect to \eqn{I_p}, both locating in Riemannian manifold \eqn{\mathcal{S}_p^+}.
#'
#' @export
#' @examples
#' C1 = CorrM(5)
#' cor.ortho_basis(C1)
cor.ortho_basis <- function(C,D = NULL){
  p = dim(C)[1]
  m = 0.5*p*(p-1)
  # optimal position to I
  if(is.null(D)) D <- D.optimal.position(diag(p),C)
  S = phi(C,D)
  # orthonormal basis in T_I of S
  #  Note that the basis of horizontal space and tangent space T_I in C are the same
  E = euc.ortho_basis(p)
  horE = E[-diag_ind(p)]
  horS = lapply(1:m, function(i) sqrtm(S)%*%horE[[i]]%*%sqrtm(S))
  # orthonormal basis in T_C
  EC = lapply(1:m, function(i) deri_submer(horS[[i]],S))
  return(EC)
}

#' Orthonormal basis in horizontal space
#' @description
#' Construct orthonormal basis in the horizontal subspace \eqn{\mathcal{H}_S}.
#'
#' @param S A SPD matrix in \eqn{\mathcal{S}_p^+}.
#'
#' @return Orthonormal basis with the dimension of \eqn{m = p(p-1)/2}
#' @export
#'
#' @examples S = CovM(5)
#' hor.ortho_basis(S)
hor.ortho_basis <- function(S){
  p = dim(S)[1]
  m = 0.5*p*(p-1)
  # orthonormal basis in T_I of S
  E = euc.ortho_basis(p)
  ### Dimension of horizontal space
  horE = E[-diag_ind(p)]
  horS = lapply(1:m, function(i) sqrtm(S)%*%horE[[i]]%*%sqrtm(S))
  return(horS)
}

#' Coordinate vector by vectorizing projection in quotient geometry
#' @description Vectorize a tangent vector matrix \eqn{X \in T_C\mathcal{S}_p^+} into a coordinate vector \eqn{v \in \mathbb{R}^m}. Thereby,
#' the dimension of affine invariant geometry \eqn{\mathcal{C}_p^+} is \eqn{m = p(p-1)/2}.
#'
#' @param X A tangent vector matrix in \eqn{T_{C}\mathcal{C}_p^+}.
#' @param C A correlation matrix in \eqn{\mathcal{C}_p^+}.
#' @param D A diagonal matrix in group action \eqn{\mathcal{D}_p^+}. It makes \eqn{C_1, DC_2D \in \mathcal{S}_p^+} being in optimal positions.
#' When it is \code{NULL}, we can obtain it from [D.optimal.position()].
#'
#' @return  A coordinate vector in \eqn{\mathbb{R}^m}
#' @export
#'
#' @examples C1 = CorrM(5) ;C2 = CorrM(5)
#' X1 = cor.LogMap(C1,C2)
#' cor.coor_vec(X1,C1)
cor.coor_vec <- function(X,C,D = NULL){
  p <- ncol(C)
  m <- 0.5*p*(p-1)
  ### optimal position w.r.t I for C
  if(is.null(D)) D <- D.optimal.position(diag(p),C)
  S = phi(C,D)
  ### horizontal lift
  horV <- horlift(X,S)
  ### orthonormal basis
  horS = hor.ortho_basis(S)
  ### orthonormal coordinate vector
  u = sapply(1:m, function(i) spd.inner_product(horV,horS[[i]],S))
  return(u)
}

#' Frechet mean by Sturm algorithm in quotient geometry
#' @description Compute Frechet sample mean in quotient geometry \eqn{\mathcal{C}_p^+}.
#'
#' @param C Array object with \eqn{p \times p \times n}.
#' @param MaxIt Maximum iterations of algorithm.
#' @param store.est Logical values with \code{FALSE} as default value. If it is \code{TRUE}, we will store the estimated means for each iteration, vice versa.
#' @param store.fsv Logical values with \code{FALSE} as default value. If it is \code{TRUE}, we will store the estimated variance for each iteration, vice versa.
#'
#' @return Frechet mean
#' @export
#'
cor.mean_Sturm = function(C,MaxIt,store.est=FALSE,store.fsv=FALSE) {
  if(is.null(store.est))  store.est = FALSE
  if(is.null(store.fsv))  store.fsv = FALSE
  # store Frechet variance
  if (store.fsv == TRUE) {
    fsv = c()
  }
  # store estimated means
  if (store.est == TRUE) {
    Sk = list()
  }
  # No. of data points
  n = dim(C)[3]
  # choose a random initial guess
  i = sample(c(1:n),1)
  old.est = C[,,i]
  cat("iterations:")
  for (k in 1:MaxIt){
    if(k == MaxIt){
      cat(k,"\n")
    } else {
      cat(k,",",sep = "")
    }
    # obtain the iterated mean based on geodesic
    j = sample(c(1:n),1)
    ranC = C[,,j]
    new.est = cor.geodesic(1/(k+1), old.est,ranC)
    if (store.est == TRUE) {
      Sk[[k]] = new.est
    }
    if (store.fsv == TRUE) {
      omega = 0
      for (i in 1:n){
        omega = omega + (cor.metric(new.est, C[,,i])^2)/n
      }
      fsv[k] = omega
    }
    old.est = new.est
  }
  if (store.est == TRUE && store.fsv == TRUE) {
    Results = list(mean = new.est, ests = Sk, frechetsv = fsv)
  } else if (store.est == TRUE && store.fsv == FALSE) {
    Results = list(mean = new.est, ests = Sk)
  } else if (store.est == FALSE && store.fsv == TRUE) {
    Results = list(mean = new.est, frechetsv = fsv)
  } else {
    Results = new.est
  }
  return(Results)
}
#'
#'
#' #' Frechet mean by Riemannian submersion
#' #' @note Parallel transport is used to obtain optimal positions for each iterations
#' #' @param C array data with \eqn{p \times p \times n}
#' #' @param MaxIt Maximum iterations of algorithm
#' #' @param store.est Logical values with FALSE as default value. If it is TRUE, we will store the estimated means for each iteration, vice versa.
#' #' @param store.fsv Logical values with FALSE as default value. If it is TRUE, we will store the estimated variance for each iteration, vice versa.
#' #' @param parallel Logical valuess with TRUE as default value. If it is TRUE, we obtain optimal positions with parallel computation
#' #' @param method.Sturm Sturm algorithm and operator from Lenglet paper are used to obtain Frechet mean in \eqn{\mathcal{S}_p^+}. Default method is "Sturm"
#' #'
#' #'
#' #' @references \insertRef{lenglet2006statistics}{geomTS}
#' #' @references \insertRef{riquelme2021riemannian}{geomTS}
#' #'
#' #' @return Frechet mean in \eqn{\mathcal{S}_p^+}
#' #' @export
#' #' @import foreach
#' #' @importFrom  doParallel registerDoParallel
#' #' @importFrom parallel detectCores
#' #'
#' cor.mean_submer <-  function(C,MaxIt, store.est=FALSE,store.fsv=FALSE,parallel = TRUE, method.Sturm = TRUE){
#'   if(is.null(store.est))  store.est = FALSE
#'   if(is.null(store.fsv))  store.fsv = FALSE
#'   if(is.null(parallel)) parallel = TRUE
#'   if(is.null(method.Sturm)) method.Sturm = TRUE
#'   # store Frechet variance
#'   if (store.fsv == TRUE) {
#'     fsv = c()
#'   }
#'   # store estimated means
#'   if (store.est == TRUE) {
#'     Sk = list()
#'   }
#'   # No. of data points
#'   n = dim(C)[3]
#'   p = dim(C)[1]
#'   # choose a random initial guess
#'   i = sample(c(1:n),1)
#'   old.est = C[,,i]
#'   ### iterations
#'   cat("iterations:")
#'   for(k in 1:MaxIt){
#'     if(k == MaxIt){
#'       cat(k,"\n")
#'     } else {
#'       cat(k,",",sep="")
#'     }
#'     if(parallel) {
#'       ### find the optimal position w.r.t C0
#'       registerDoParallel(detectCores())
#'       newS <- foreach(j = 1:n,.packages = "expm")%dopar%{
#'         # cat("j = ",j,"\n")
#'         D <- D.optimal.position(old.est,C[,,j])
#'         tmp <- phi(C[,,j],D)
#'         return(tmp)
#'       }
#'       newS = list2array(newS)
#'     } else {
#'       newS  = array(0,c(p,p,n))
#'       for(j in 1:n){
#'         # cat("j = ",j,"\n")
#'         D <- D.optimal.position(old.est,C[,,j])
#'         newS[,,j] <- phi(C[,,j],D)
#'       }
#'     }
#'     if(method.Sturm){
#'       ### compute the mean for newS by Sturm algorithm
#'       FM_S <- spd.mean_Sturm(newS,MaxIt = 200)
#'     } else {
#'       ### compute the mean for newS by operator from Lenglet paper
#'       FM_S <- spd.mean_Lenglet(newS,MaxIt = 200,dt = 0.05)
#'     }
#'     ### project back to CORR
#'     new.est <- submer(FM_S)
#'     if (store.est == TRUE) {
#'       Sk[[k]] = new.est
#'     }
#'     if (store.fsv == TRUE) {
#'       omega = 0
#'       for (i in 1:n){
#'         omega = omega + (cor.metric(new.est, C[,,i])^2)/n
#'       }
#'       fsv[k] = omega
#'     }
#'     old.est = new.est
#'   }
#'   if (store.est == TRUE && store.fsv == TRUE) {
#'     Results = list(mean = new.est, ests = Sk, frechetsv = fsv)
#'   } else if (store.est == TRUE && store.fsv == FALSE) {
#'     Results = list(mean = new.est, ests = Sk)
#'   } else if (store.est == FALSE && store.fsv == TRUE) {
#'     Results = list(mean = new.est, frechetsv = fsv)
#'   } else {
#'     Results = new.est
#'   }
#'   return(Results)
#' }
