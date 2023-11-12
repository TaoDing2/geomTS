#' Simulate SPD matrices under the framework of manifold-adapted time series model with scalar coefficients
#' @description
#'  Simulating SPD matrices in affine invariant manifold \eqn{\mathcal{S}_p^+} according to the mean-reverting autoregressive model
#' with scalar coefficients.
#'
#' @param n Number of total simulated time points.
#' @param alpha A vector. It represents the autoregressive coefficients.
#' @param sigma Standard deviation of white noise
#' @param beta The regression coefficient of mean reverting term. It should be given when \code{attS} is not \code{NULL}.
#' @param attS An attract point in \eqn{\mathcal{S}_p^+}. It should be given when \code{beta} is not \code{NULL}.
#' @param p The dimension \eqn{p}, if \code{beta} and \code{attS} are missing.
#'
#' @return An  \eqn{p \times p \times n} array data in \eqn{\mathcal{S}_p^+}.
#' @export
#'
spd.simulate <- function(n,alpha,sigma,beta = NULL,attS= NULL,p = NULL){
  if(is.null(p)){
    if(!is.null(attS)) {
      p = dim(attS)[1]
    } else {
      stop("Please input the dimension of p")
    }
  }
  m   <- 0.5*p*(p+1)  # dimensions
  q = length(alpha)   # lags
  ###
  ### initial points
  ###
  # We generate the initial data points from [CovM()]
  S <- array(0,c(p,p,n))
  for(i in 1:(q+1)){
    S[,,i] <- CovM(p)
  }
  ### tangent vectors
  newV <- list()
  for(i in 1:q){
    newV[[i]] <- spd.LogMap(S[,,i],S[,,i+1])
  }
  ### coordinate vectors
  newv <- matrix(0, ncol= m, nrow = n-1)
  for(i in 1:q){
    newv[i,] <- spd.coor_vec(newV[[i]],S[,,i])
  }
  # orthonormal basis on Tangent space T_{I}
  E <- euc.ortho_basis(p)
  ###
  ### Simulate data
  ###
  for(i in (q+1):(n-1)){
    if(i == q+1){
      cat("Simulating i = ",i,",",sep = "")
    } else if(i == (n-1)){
      cat(i,"\n")
    } else {
      cat(i,",",sep = "")
    }
    mu <- rep(0,length=m)
    # sum of multiplication of lag vectors and regression regression
    for(l in 1:q){
      Vtmp <- spd.para_tran(S[,,i-l],S[,,i],newV[[i-l]])
      mu <- mu + alpha[l]*spd.coor_vec(Vtmp,S[,,i])
    }
    # mean reversion
    if(!is.null(beta)){
      mu = mu + beta*spd.coor_vec(spd.LogMap(S[,,i],attS),S[,,i])
    }
    # Euclidean vector
    newv[i,] <- rnorm(m,mu,sigma)
    # project back to tangent space
    matV <- matrix(0,p,p)
    for(r in 1:m){
      matV <- matV + newv[i,r] * sqrtm(S[,,i])%*% E[[r]] %*%sqrtm(S[,,i])
    }
    newV[[i]] <- matV
    # Exponential map
    S[,,i+1] = spd.ExpMap(S[,,i],matV)
  }
  return(S)
}

#' Simulate symmetric matrices in \eqn{Sym(p)}
#'
#' @param n Number of time points
#' @param alpha Coefficients in VAR model whose length is the lags of VAR terms
#' @param sigma Standard deviation of white noise
#' @param beta Coefficient of mean reversion. It follows the attract point/ Mean
#' @param attS Attract point.
#' @param p Dimensions if \eqn{beta} and \eqn{attS} are missing.

#'
#' @return A array data \eqn{p \times p \times n}
#' @export
#'
euc.simulate <- function(n,alpha,sigma,beta = NULL,attS = NULL,p = NULL){
  if(is.null(attS)) attS = diag(p)
  if(is.null(p)){
    p = dim(attS)[1]
  }
  m    <- 0.5*p*(p+1)     # dimensions
  q = length(alpha)
  ###
  ### initial points
  ###
  S <- array(0,c(p,p,n))
  for(i in 1:(q+1)){
    S[,,i] <- CovM(p)
  }
  ### Initial tangent vectors
  newV <- list()
  for(i in 1:q){
    newV[[i]] <- euc.LogMap(S[,,i],S[,,i+1])
  }
  ### coordinate vectors
  newv <- matrix(0, ncol= m, nrow = n-1)
  for(i in 1:q){
    newv[i,] <- euc.coor_vec(newV[[i]])
  }
  # orthonormal basis on Tangent space T_{I}
  E <- euc.ortho_basis(p)
  ###
  ### Simulate data
  ###
  for(i in (q+1):(n-1)){
    if(i == q+1){
      cat("Simulating i = ",i,",",sep = "")
    } else if(i == (n-1)){
      cat(i,"\n")
    } else {
      cat(i,",",sep = "")
    }
    # sum of multiplication of lag vectors and regression regression
    mu <- rep(0,length=m)
    for(l in 1:q){
      mu <- mu + alpha[l]*euc.coor_vec(newV[[i-l]])
    }
    # mean reversion
    if(!is.null(beta)){
      mu = mu + beta*euc.coor_vec(euc.LogMap(S[,,i],attS))
    }
    # Euclidean vector
    newv[i,] <- rnorm(m,mu,sigma)
    # project back to tangent space
    matV <- matrix(0,p,p)
    for(r in 1:m){
      matV <- matV + newv[i,r] *E[[r]]
    }
    newV[[i]] <- matV
    # Exponential map
    S[,,i+1] <- euc.ExpMap(S[,,i],matV)
  }
  return(S)
}
