################################################################################
###
###  Affine invariant manifold consisting of symmetric positive definite (SPD) matrices
###
################################################################################
#' Geodesic distance in affine invariant geometry
#' @description Calculate geodesic distance for any two SPD matrices \eqn{S_1,S_2} in SPD matrix space \eqn{\mathcal{S}_p^+}.
#'
#' @param S1 A SPD matrix.
#' @param S2 A SPD matrix.
#'
#' @return A numerical value of geodesic distance between \eqn{S_1} and \eqn{S_2} in \eqn{\mathcal{S}_p^+}
#' @export
#'
#' @details
#' Given any two points \eqn{S_1} and \eqn{S_2}  in SPD matrix space \eqn{\mathcal{S}_p^+}, we can define the Riemannian affine invariant metric as
#' \deqn{d(S_1,S_2) = \sqrt{\sum_{i=1}^{p}\log^2(\lambda_i)}}
#' where \eqn{\lambda_i,i = 1,\ldots,p} are eigenvalues of matrix \eqn{S_1^{-1/2}S_2S_1^{-1/2}}.
#'
#' @note Note that \eqn{d(S_1,S_2) = d(S_2,S_1)}
#' @examples S1 = CovM(3)
#' S2 = CovM(3)
#' spd.metric(S1,S2)
#' spd.metric(S2,S1)
spd.metric <- function(S1,S2){
  A <- solve(sqrtm(S1))
  tmp <- A%*%S2%*%A
  eta <- eigen(tmp)$values
  dist <- sqrt(sum(log(eta)^2))
  return(dist)
}

#' Inner product in affine invariant geometry
#' @description Compute the inner product for any two tangent vector matrices \eqn{V_1, V_2 \in T_S\mathcal{S}_p^+} at the base point \eqn{S \in \mathcal{S}_p^+ }.
#'
#' @param V1 A tangent vector matrix in  \eqn{T_S\mathcal{S}_p^+}.
#' @param V2 A tangent vector matrix in  \eqn{T_S\mathcal{S}_p^+}.
#' @param S A SPD matrix. It has a tangent space \eqn{T_S\mathcal{S}_p^+}.
#'
#' @return A numerical distance.
#' @export
#'
#' @details For any tangent vectors \eqn{V_1, V_2} in the tangent space \eqn{T_S\mathcal{S}_p^+} at the base point of \eqn{S \in \mathcal{S}_p^+},
#' the inner product of the Riemannian metric is defined as
#' \deqn{ \langle V_1,W_2 \rangle_S =  \mathrm{tr}(S^{-1}V_1S^{-1}W_2).}
#'
#' @note Two tangent vectors should both lie in the same tangent space \eqn{T_S\mathcal{S}_p^+}.
#' @examples S1 = CovM(3)
#' S2 = CovM(3)
#' S3 = CovM(3)
#' V1 = spd.LogMap(S1,S2)
#' V2 = spd.LogMap(S1,S3)
#' spd.inner_product(V1,V2,S1)
spd.inner_product <- function(V1,V2,S){
  invs <- solve(S)
  return(tr(invs%*%V1%*%invs%*%V2))
}

#' Logarithm map in affine invariant geometry
#' @description Logarithm map in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#' It aims to obtain a tangent vector matrix by computing the velocity from point \eqn{S_1} to point \eqn{S_2} both locating in \eqn{\mathcal{S}_p^+}.
#'
#' @param S1 A SPD matrix.
#' @param S2 A SPD matrix.
#'
#' @return A tangent vector matrix in \eqn{T_{S}\mathcal{S}_p^+}.
#'
#' @note \itemize{
#' \item \eqn{T_{S}\mathcal{S}_p^+} is equivalent to the symmetric matrix space \eqn{Sym(p)}.
#' \item The tangent vector is the initial velocity at the starting point in the geodesic connecting \eqn{S_1} and \eqn{S_2} both in\eqn{\mathcal{S}_p^+}.
#' \item It is the inverse map of [spd.ExpMap()].
#' }
#'
#' @export
#' @seealso [spd.ExpMap()]
#' @examples S1 = CovM(3);S2 = CovM(3)
#' spd.LogMap(S1,S2)
spd.LogMap <- function(S1,S2){
  sqt  <- sqrtm(S1)
  invS <- solve(sqt)
  V <- sqt%*%logm(invS%*%S2%*%invS)%*%sqt
  return(V)
}

#' Exponential map in affine invariant geometry
#' @description Exponential map in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#' It projects the tangent vector matrix \eqn{V \in T_{S}\mathcal{S}_p^+} at point \eqn{S \in \mathcal{S}_p^+} back to \eqn{\mathcal{S}_p^+}.
#'
#' @param S A SPD matrix.
#' @param V A tangent vector matrix in \eqn{T_{S}\mathcal{S}_p^+}.
#'
#' @return A SPD matrix. It is the mapped point of \eqn{V} at base point \eqn{S} in the manifold \eqn{\mathcal{S}_p^+}.
#' @export
#'
#' @note It is the inverse map of [spd.LogMap()].
#'
#' @seealso [spd.LogMap()]
#' @examples S1 = CovM(3)
#' S2 = CovM(3)
#' V1 = spd.LogMap(S1,S2)
#' spd.ExpMap(S1,V1)
spd.ExpMap <- function(S,V){
  sqt  <- sqrtm(S)
  invs <- solve(sqt)
  return(sqt%*%expm(invs%*%V%*%invs)%*%sqt)
}


#' @title Geodesic in affine invariant geometry
#' @description Compute geodesic \eqn{\gamma(t), t \in [0,1]} in affini invariant geometry \eqn{\mathcal{S}_p^+}.
#'
#' @param t Time. It should be in the interval of \eqn{[0,1]}.
#' @param S1 A SPD matrix. It is the starting point of the geodesic  \eqn{\gamma(t =0)} in \eqn{\mathcal{S}_p^+}.
#' @param S2 A SPD matrix. It is the endpoint of the geodesic \eqn{\gamma(t =1)} in \eqn{\mathcal{S}_p^+}.
#'
#' @details It is defined by \eqn{\mathrm{Exp}_{S_1}(t V)}, where \eqn{V = \mathrm{Log}_{S_1}(S_2)} and \eqn{t \in [0,1]}.
#'
#' @return A SPD matrix. It lies on the geodesic \eqn{\gamma(t)} at time \eqn{t, t\in[0,1]}.
#' @export
#'
#' @examples S1= CovM(3)
#' S2 = CovM(3)
#' spd.geodesic(t = 0.5,S1,S2)
#' spd.geodesic(t = 0.1,S1,S2)
spd.geodesic <- function(t,S1,S2){
  V1 = spd.LogMap(S1,S2)
  sqt = sqrtm(S1)
  invs = solve(sqt)
  geodesic = sqt%*%expm(t * invs%*%V1%*%invs)%*%sqt
  return(geodesic)
}


#' @title Parallel transport in affine invariant geometry
#' @description Parallel transport the tangent vector \eqn{V \in T_{S_1}\mathcal{S}_p^+}
#' into another tangent space \eqn{T_{S_2}\mathcal{S}_p^+} along the geodesic connecting \eqn{S_1} and \eqn{S_2} in \eqn{\mathcal{S}_p^+}.
#'
#' @param S1 A SPD matrix.
#' @param S2 A SPD matrix.
#' @param V A tangent vector matrix in \eqn{T_{S_1}\mathcal{S}_p^+}.
#'
#' @details  Parallel transport involves the translation of a tangent vector \eqn{V} from the tangent space \eqn{T_{S_1}\mathcal{S}_p^+}
#'  at point \eqn{S_1 \in \mathcal{S}_p^+} to the tangent space \eqn{T_{S_2}\mathcal{S}_p^+} at point \eqn{S_2 \in \mathcal{S}_p^+}
#'  and its  mathematical formula is :
#'  \deqn{\mathcal{P}_{S_1\rightarrow S_2}(V) = WV_1W^T.}
#' Here, \eqn{W} can be simply written  as  \eqn{W = (S_2S_1^{-1})^{1/2}}.
#'
#' @references \insertRef{yair2019parallel}{geomTS}
#'
#' @return A translated tangent vector matrix in \eqn{T_{S_2}\mathcal{S}_p^+}.
#' @export
#'
#' @examples S1 = CovM(3)
#' S2 = CovM(3)
#' V1 = spd.LogMap(S1,S2)
#' spd.para_tran(S1,S2,V1)
spd.para_tran <- function(S1,S2,V){
  W <- sqrtm(S2%*%solve(S1))
  PT <-  W%*%V%*%t(W)
  return(PT)
}


#' @title Orthonormal basis at any base point in affine invariant geometry
#' @description Construct orthonormal basis in the tangent space \eqn{T_S\mathcal{S}_p^+}.
#'
#' @param S A SPD matrix.
#'
#' @return A list of orthonormal basis in the tangent space \eqn{T_S\mathcal{S}_p^+}
#' containing \eqn{m = \frac{1}{2}p(p+1)} observations of \eqn{p \times p} basis.
#'
#' @details An orthonormal basis of \eqn{T_S\mathcal{S}_p^+} for each  \eqn{S \in \mathcal{S}_p^+} can be defined
#' via parallel transport of a choice of orthonormal basis of the tangent space at the identity matrix \eqn{I_p \in T_{I_p}\mathcal{S}_p^+} .
#' The matrices \eqn{E_r, r= 1,\ldots,m = p(p+1)/2}, defined in [euc.coor_vec()] determine an orthonormal basis of the tangent space at
#' \eqn{I_p} with respect to the Riemannian inner product. Using parallel transport from \eqn{I_p} to a general point \eqn{S \in \mathcal{S}_p^+} in [spd.para_tran()],
#' we obtain an orthonormal basis of \eqn{T_S\mathcal{S}_p^+} defined by
#' \deqn{ E^{S}_r = \mathcal{P}_{I_p \rightarrow S}(E_r) = S^{1/2}E_rS^{1/2}. }
#' This basis can be thought of as a global smooth section of the frame bundle of the Riemannian manifold \eqn{\mathcal{S}_p^+}.
#'
#' @note By definition, it follows that affine invariant geometry is parallelizable.

#' @export
#'
#' @examples
#' spd.ortho_basis(diag(3))
#' S1 = CovM(3)
#' spd.ortho_basis(S1)
spd.ortho_basis <- function(S){
  p <- ncol(S)
  m = 0.5*p*(p+1)
  ###
  ### orthonormal basis on T_{I}
  ###
  E <- euc.ortho_basis(p)
  ###
  ###  Return the orthonormal basis on T_{S}
  ###
  ES <- list()
  sqt <- sqrtm(S)
  for(i in 1:m){
    ES[[i]] <-  sqt%*%E[[i]]%*%sqt
  }
  return(ES)
}


#' Coordinate vector by vectorizing projection in affine invariant geometry
#' @description Vectorize a tangent vector matrix \eqn{V \in T_S\mathcal{S}_p^+} into a coordinate vector \eqn{v \in \mathbb{R}^m}. Thereby,
#' the dimension of affine invariant geometry \eqn{\mathcal{S}_p^+} is \eqn{m = p(p+1)/2}.
#'
#' @param V A tangent vector matrix in \eqn{T_S\mathcal{S}_p^+} .
#' @param S A SPD matrix.
#'
#' @details Given any tangent vector \eqn{V \in T_S\mathcal{S}_p^+},
#' we can calculate its orthonormal coordinate vector \eqn{\boldsymbol{v} = (v_1,...,v_m)^T} in \eqn{\mathbb{R}^m} as
#' \deqn{\boldsymbol{v} = \langle V,S^{1/2}\boldsymbol{E}S^{1/2}  \rangle_S = \mathrm{tr}( S^{-1/2}VS^{-1/2}\boldsymbol{E}))}
#'  where \eqn{\boldsymbol{E} = (E_1,\ldots,E_m)} is the orthonormal basis in \eqn{I_p} as shown in [spd.ortho_basis].
#'  Thereby, vectorization of \eqn{V \in T_S\mathcal{S}_p^+}  is mapped to \eqn{\boldsymbol{v} = (v_1,\ldots,v_m)^T \in \mathbb{R}^m}.
#'
#' @return A coordinate vector in \eqn{\mathbb{R}^m}.
#' @export
#' @seealso [spd.Euc2Manif()]
#' @examples S1 = CovM(3)
#' S2 = CovM(3)
#' V1 = spd.LogMap(S1,S2)
#' spd.coor_vec(V1,S1)
spd.coor_vec <- function(V,S){
  p <- ncol(S)
  m <- 0.5*p*(p+1)
  S_basis <- spd.ortho_basis(S) # orthonormal basis on T_{S}
  v <- sapply(1:m, function(i) spd.inner_product(V,S_basis[[i]],S))
  return(v)
}

#' Mapping back the coordinate vector to data point in affine invariant geometry
#' @description
#' Map back the coordinate vector \eqn{\boldsymbol{v} \in \mathbb{R}^m} to the tangent space \eqn{T_S\mathcal{S}_p^+}
#' and then do exponential map the base point \eqn{S \in \mathcal{S}_p^+}.
#'
#' @param v Coordinate vector in \eqn{\mathbb{R}^m}.
#' @param S A SPD matrix.
#'
#' @details It follows the following process:
#' \deqn{V = \sum_{r = 1}^{m}v_rS^{1/2}E_rS^{1/2}\; \text{and then} \; S^* = \mathrm{Exp}_{S}(V).}
#'
#' @return A SPD matrix in \eqn{\mathcal{S}_p^+}
#' @export
#' @seealso [spd.coor_vec()]
#' @examples
#' S1 = CovM(5); S2 = CovM(5)
#' V = spd.LogMap(S1,S2)
#' v = spd.coor_vec(V,S1)
#' spd.Euc2Manif(v,S1)
spd.Euc2Manif <- function(v,S){
  # dimension
  m = length(v)
  p = ncol(S)
  # basis in tangent space
  E = spd.ortho_basis(S)
  ### Tangent vector
  V = matrix(0,p,p)
  for(i in 1:m){
    V = V + v[i]*E[[i]]
  }
  ### point in manifold
  return(spd.ExpMap(S,V))
}

#' Frechet mean in affine invariant geometry
#' @description
#' Compute Frechet sample mean in affine invariant geometry \eqn{\mathcal{S}_p^+}
#'
#' @param S Array object with \eqn{p \times p \times n}.
#' @param MaxIt Maximum iterations of algorithm
#' @param conver Threshold of convergence. Default value is \eqn{10^{-4}}
#' @param dt Step size of algorithm. Normally, the default value is set  as 0.1.
#' @param store.M Logical values with \code{FALSE} as default value. If it is \code{TRUE}, we will store the estimated means and variance for each iteration, vice versa.
#' @param method Either "Sturm" (Sturm algorithm) nor "GD" (Graident descent algorithm).
#' @references \insertRef{lenglet2006statistics}{geomTS}
#' \insertRef{sturm2003probability}{geomTS}
#'
#' @return Estimated Frechet mean (and Frechet variance at each iteration if store.M = \code{TRUE})
#' @export
#'
#' @examples
#' S = lapply(1:10, function(i) CovM(5))
#' S = list2array(S)
#' FM = spd.mean(S,method = "Sturm",store.M = TRUE)
#' FM$mean
#'
spd.mean = function(S,MaxIt = NULL,conver = NULL,dt = NULL,store.M = NULL, method = c("Sturm", "GD")) {
  if(is.null(store.M))  store.M = FALSE
  if(is.null(MaxIt)) MaxIt = 200
  if(is.null(conver)) conver = 4
  if(is.null(dt))  dt = 0.1
  ### store estimated means
  if (store.M == TRUE) {
    itS = list()
  }
  ### Objective function
  fsv <- c()
  ### No. of data points
  n = dim(S)[3]
  ### The initial guess
  S0 = S[,,sample(1:n,1)]
  ###
  ### iterations
  ###
  k = 1
  convergence = 1
  while (convergence > 10^{-conver}) {
    if(method == "GD"){
      ### compute the operator
      p = dim(S)[1]
      sqt = sqrtm(S0)
      invS = solve(sqt)
      sumM = matrix(0,p,p)
      for (j in 1:n){
        sumM = sumM + logm(invS %*% S[,,j] %*% invS)
      }
      S1 = sqt%*% expm((dt/n) * sumM)%*%sqt
    } else if ( method == "Sturm") {
      ### A random point in data set
      ranS = S[,,sample(1:n,1)]
      S1 = spd.geodesic(1/(k+1), S0,ranS)
    } else {
      print("PLease choose a method")
      break
    }
    ### store updated means if required
    if (store.M == TRUE) {
      itS[[k]] = S1
    }
    ### Frechet variance
    fsv[k] = mean(sapply(1:n, function(j) spd.metric(S1,S[,,j])^2))
    ### iterating criterion
    if(k > 1){
      convergence = abs(fsv[k] - fsv[k-1])
    }
    ### forwarding proceeds
    k = k+1
    S0 = S1
    if(k > MaxIt) break
  }
  if (store.M) {
    results = list(mean = S0, itS = itS, fsv = fsv)
  } else {
    results = S0
  }
  return(results)
}

