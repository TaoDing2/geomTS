#' Exponential map in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#' @description
#' Exponential map from \eqn{S} into direction \eqn{V}
#'
#'
#' @param S Base point in manifold \eqn{\mathcal{S}_p^+}
#' @param V Tangent vector in tangent space \eqn{T_S\mathcal{S}_p^+}
#'
#' @return Mapped point in \eqn{\mathcal{S}_p^+}
#' @export
#'
lch.ExpMap <- function(S,V){
  x = t(chol(S))
  tmp = solve(x)%*%V%*%t(solve(x))
  ### extract lower triangular matrix with diagonal multiplied by 0.5
  tmp[upper.tri(tmp,diag = FALSE)] = 0
  W = x%*%(tmp - 0.5* Diag(tmp))
  ### exponential map on Symmetric Space
  tmp = cho.Exp(x,W)
  return(tmp%*%t(tmp))
}


#' Logarithm map in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#' @description
#' Logarithm map emanating from \eqn{S_1} to \eqn{S_2}
#'
#' @param S1 Starting point in manifold \eqn{\mathcal{S}_p^+}
#' @param S2 Ending point in manifold \eqn{\mathcal{S}_p^+}
#'
#' @return Tangent vector in \eqn{T_S\mathcal{S}_p^+}
#' @export
#'
lch.LogMap <- function(S1,S2){
  x = t(chol(S1))
  y = t(chol(S2))
  W = cho.Log(x,y)
  tmp = x%*%t(W) + W%*%t(x)
  return(tmp)
}

#' Geodesic in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#'
#' @param t Time which should be in the interval of \eqn{[0,1]}
#' @param S1 Starting point in manifold \eqn{\mathcal{S}_p^+}
#' @param S2 Ending point in manifold \eqn{\mathcal{S}_p^+}
#'
#' @return Data ponint the geodesic from \eqn{S_1} to \eqn{S_2}
#' @export
#'
lch.geodesic <- function(t,S1,S2){
  X = lch.LogMap(S1,S2)
  return(lch.ExpMap(S1,t*X))
}



#' Geodesic distance in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#'
#' @param S1 Starting point in manifold \eqn{\mathcal{S}_p^+}
#' @param S2 Ending point in manifold \eqn{\mathcal{S}_p^+}
#'
#' @return A numerical value of distance between \eqn{S_1} and \eqn{S_2} in \eqn{\mathcal{S}_p^+}
#' @export
#'
lch.metric <- function(S1,S2){
  x = t(chol(S1))
  y = t(chol(S2))
  return(cho.metric(x,y))
}

#' Inner product in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#'
#' @param X Tangent vector in tangent space \eqn{T_S\mathcal{S}_p^+}
#' @param Y Tangent vector in tangent space \eqn{T_S\mathcal{S}_p^+}
#' @param S Base point
#'
#' @return A numerical value which measures the inner product between two tangent vectors
#' @export
#'
lch.inner_product <- function(X,Y,S){
  z = t(chol(S))
  ### extract lower triangular matrix with diagonal multiplied by 0.5
  tmp = solve(z)%*%X%*%t(solve(z))
  tmp[upper.tri(tmp,diag = FALSE)] = 0
  W1 = z%*%(tmp - 0.5* Diag(tmp))
  ###
  tmp = solve(z)%*%Y%*%t(solve(z))
  tmp[upper.tri(tmp,diag = FALSE)] = 0
  W2 = z%*%(tmp - 0.5* Diag(tmp))
  ###
  return(cho.inner_product(W1,W2,z))
}


#' Parallel transport in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#'
#' @param S1 Base point in \eqn{\mathcal{S}_p^+}
#' @param S2 Objective point in \eqn{\mathcal{S}_p^+}
#' @param X Tangent vector in tangent space \eqn{T_{S_1}\mathcal{S}_p^+}
#'
#' @return Transported vector along the geodesic at \eqn{T_{S_2}\mathcal{S}_p^+}
#' @export
#'
lch.para_tran <- function(S1,S2,X){
  x = t(chol(S1))
  y = t(chol(S2))
  ### extract lower triangular matrix with diagonal multiplied by 0.5
  tmp = solve(x)%*%X%*%t(solve(x))
  tmp[upper.tri(tmp,diag = FALSE)] = 0
  W = x%*%(tmp - 0.5* Diag(tmp))
  V = cho.para_tran(x,y,W)
  return(y%*%t(V) + V%*%t(y))
}


#' Parallel transport in \eqn{\mathcal{S}_p^+} with LogCholesky metric
#'
#' @param S Array data
#'
#' @return Mean
#' @export
#'
lch.mean <- function(S){
  n = dim(S)[3]
  L = lapply(1:n, function(i) t(chol(S[,,i])))
  M = cho.mean(list2array(L))
  return(M%*%t(M))
}





