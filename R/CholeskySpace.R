#' Exponential map on the Cholesky Space
#' @description
#' Compute expoenential map on the Cholesky Space emanating from the lower triangular
#' matrix with positive diagonal \eqn{L} towards the lower triangular matrix \eqn{X}
#'
#'
#' @param L lower triangular matrix with positive diagonal
#' @param X lower triangular matrix
#'
#' @return Lower triangular matrix after exponential map
#' @export
#'
cho.Exp <- function(L,X){
  tmp = L - Diag(L) + X - Diag(X)  + Diag(L)%*%expm(Diag(X)%*%solve(Diag(L)))
  return(tmp)
}

#' Logarithmic map on the CholeskySpace
#' @description
#' Compute the logarithmic map on the Cholesky Space for  the geodesic emanating from
#'  the lower triangular matrix with positive diagonal \eqn{L_1} towards \eqn{L_2}.
#'
#' @param L1 lower triangular matrix with positive diagonal
#' @param L2 lower triangular matrix with positive diagonal
#'
#' @return Lower triangular matrix after logarithm map
#' @export
#'
cho.Log <- function(L1,L2){
  tmp = L2 - L1 - Diag(L2) + Diag(L1)  + Diag(L1)%*%logm(solve(Diag(L1))%*%Diag(L2))
  return(tmp)
}

#' Geodesic on the Cholesky Space
#' @description
#' Geodesic emanating from p with direction X
#'
#'
#' @param t time \eqn{t, t \in [0,1]}
#' @param L1 lower triangular matrix
#' @param L2 lower triangular matrix
#'
#' @return Geodesic
#' @export
cho.geodesic <- function(t,L1,L2){
  X = cho.Log(L1,L2)
  tmp = L1 - Diag(L1) + t(X - Diag(X))  + Diag(L1)%*%expm(t*Diag(X)%*%solve(Diag(L1)))
  return(tmp)
}

#' Geodesic distance on the Cholesky Space between two matrices \eqn{L_1}, \eqn{L_2}
#' @description
#' Riemannian distance on the Cholesky Space between two matrices \eqn{L_1} and \eqn{L_2}
#' that are lower triangular with positive diagonal
#'
#'
#' @param L1 Lower triangular with positive diagonal
#' @param L2 Lower triangular with positive diagonal
#'
#' @return Riemannian distance
#' @export
#'
cho.metric <- function(L1,L2){
  tmp1 = norm(L1-L2-Diag(L1)+Diag(L2),type="F")^2
  tmp2 = norm(logm(Diag(L1))-logm(Diag(L2)),type= "F")^2
  return(sqrt(tmp1 + tmp2))
}


lch.metric <- function(S1,S2){
  x = t(chol(S1))
  y = t(chol(S2))
  tmp1 = norm(x-y-Diag(x)+Diag(y),type="F")^2
  tmp2 = norm(logm(Diag(x))-logm(Diag(y)),type="F")^2
  return(sqrt(tmp1+tmp2))
}


#' Inner product on the Cholesky Space
#'
#' @param X Tangent vector which is lower triangular matrix  with arbitrary diagonal
#' @param Y Tangent vector which is lower triangular matrix  with arbitrary diagonal
#' @param L Lower triangular matrix with positive diagonal
#'
#' @return Inner product
#' @export
#'
cho.inner_product <- function(X,Y,L){
  tmp1 = sum((X - Diag(X))*(Y - Diag(Y)))
  tmp2 = sum(diag(X)*diag(Y)/(diag(L)^2))
  return(tmp1 + tmp2)
}


#' Parallel transport on the Cholesky Space
#' @description
#' Parallel transport the tangent vector \eqn{X} at \eqn{L_1} along the geodesic to \eqn{L_2} on
#' Cholesky Space
#'
#' @param L1 Lower triangular with positive diagonal
#' @param L2 Lower triangular with positive diagonal
#' @param X Tangent vector which is lower triangular matrix  with arbitrary diagonal
#'
#' @return Transported vector \eqn{X} along the geodesic
#' @export
cho.para_tran <- function(L1,L2,X){
  tmp = X - Diag(X) + Diag(L2)%*%solve(Diag(L1))%*%Diag(X)
  return(tmp)
}

#' Mean on the Cholesky Space
#'
#' @param L Array data
#'
#' @return Mean
#' @export
#'
cho.mean <- function(L){
  n = dim(L)[3]
  L0  = lapply(1:n, function(i) L[,,i] - Diag(L[,,i]))
  tmp1 = euc.mean(list2array(L0))
  L1 = lapply(1:n, function(i) logm(Diag(L[,,i])))
  tmp2 = euc.mean(list2array(L1))
  return(tmp1 + expm(tmp2))
}


#' Orthonormal basis at \eqn{I_p} on the Cholesky Space
#'
#' @param p dimensions
#'
#' @return Orthomormal basis and its dimension is \eqn{m = p(p+1)/2}
#' @export
#'
cho.ortho_basis <- function(p){
  E <- list()
  k = 1
  for(i in 1:p){
    for(j in 1:i){
      ei  = ej = rep(0,p)
      ei[i] = 1
      ej[j] = 1
      E[[k]] = ei%*%t(ej)
      k = k+1
    }
  }
  return(E)
}
