#' Simulate symmetric matrices under the framework of manifold-adapted model
#' @description
#' Simulate symmetric matrices under framework of manifold-adapted model in \eqn{Sym(p)}.
#' @param n Number of simulated matrix-valued time points.
#' @param alpha Coefficients in the model whose length is the lags of autoregressive terms.
#' If it is a vector, the manifold-adapted model with scalar coefficients is used.
#' If it is a list with multiple diagonals, the model with diagonal coefficients is used.
#' @param sigma Standard deviation of white noise.
#' @param beta Coefficient of mean reversion. It follows the attract point/ Mean.
#' When we consider the mean-reverting term in the model, \code{attS} should be input and then \code{beta} is required.
#' @param attS Attract point.
#' @param p Dimensions if \eqn{beta} and \eqn{attS} are missing.

#'
#' @return An array data \eqn{p \times p \times n} in  \eqn{Sym(p)}.
#' @export
#' @importFrom MASS mvrnorm
#'
euc.simulate <- function(n,alpha,sigma,beta = NULL,attS = NULL,p = NULL){
  if(is.null(attS)) attS = diag(p)
  if(is.null(p)){
    p = dim(attS)[1]
  }
  m <- 0.5*p*(p+1)     # dimensions
  q = length(alpha)
  if(is.list(alpha)) {
    Sca = FALSE
  } else {
    Sca = TRUE
  }
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
      if(Sca){
        mu <- mu + alpha[l]*euc.coor_vec(newV[[i-l]])
      } else{
        mu <- mu + alpha[[l]]*euc.coor_vec(newV[[i-l]])
      }
    }
    # mean reversion
    if(!is.null(beta)){
      mu = mu + beta*euc.coor_vec(euc.LogMap(S[,,i],attS))
    }
    # Euclidean vector
    if(Sca){
      newv[i,] <- rnorm(m,mu,sigma)
    } else {
      newv[i,] <- mvrnorm(1,mu,diag(sigma^2))
    }
    # project coordinate vector back to manifold
    S[,,i+1] = euc.Euc2Sym(newv[i,],S[,,i],corMat = FALSE)
  }
  return(S)
}

#' Simulate SPD matrices under the framework of manifold-adapted model
#' @description
#' Simulate SPD matrices under framework of manifold-adapted model in \eqn{\mathcal{S}_p^+}.
#' @param n Number of simulated matrix-valued time points.
#' @param alpha Coefficients in the model whose length is the lags of autoregressive terms.
#' If it is a vector, the manifold-adapted model with scalar coefficients is used.
#' If it is a list with multiple diagonals, the model with diagonal coefficients is used.
#' @param sigma Standard deviation of white noise.
#' @param beta Coefficient of mean reversion. It follows the attract point/ Mean.
#' When we consider the mean-reverting term in the model, \code{attS} should be input and then \code{beta} is required.
#' @param attS Attract point.
#' @param p Dimensions if \eqn{beta} and \eqn{attS} are missing.

#'
#' @return An array data \eqn{p \times p \times n} in  \eqn{\mathcal{S}_p^+}.
#' @export
#' @importFrom MASS mvrnorm
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
  if(is.list(alpha)) {
    Sca = FALSE
  } else {
    Sca = TRUE
  }
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
      if(Sca){
        mu <- mu + alpha[l]*spd.coor_vec(Vtmp,S[,,i])
      } else{
        mu <- mu + alpha[[l]]*spd.coor_vec(Vtmp,S[,,i])
      }
    }
    # mean reversion
    if(!is.null(beta)){
      mu = mu + beta*spd.coor_vec(spd.LogMap(S[,,i],attS),S[,,i])
    }
    # Euclidean vector
    if(Sca){
      newv[i,] <- rnorm(m,mu,sigma)
    } else {
      newv[i,] <- mvrnorm(1,mu,diag(sigma^2))
    }
    # project coordinate vector back to manifold
    S[,,i+1] = spd.Euc2Manif(newv[i,],S[,,i])
  }
  return(S)
}


#' Simulate full rank correlation matrices under the framework of manifold-adapted model
#' @description
#' Simulate full rank correlation matrices under framework of manifold-adapted model in \eqn{\mathcal{C}_p^+}.
#' @param n Number of simulated matrix-valued time points.
#' @param alpha Coefficients in the model whose length is the lags of autoregressive terms.
#' If it is a vector, the manifold-adapted model with scalar coefficients is used.
#' If it is a list with multiple diagonals, the model with diagonal coefficients is used.
#' @param sigma Standard deviation of white noise.
#' @param beta Coefficient of mean reversion. It follows the attract point/ Mean.
#' When we consider the mean-reverting term in the model, \code{attC} should be input and then \code{beta} is required.
#' @param attC Attract point.
#' @param p Dimensions if \eqn{beta} and \eqn{attC} are missing.

#'
#' @return An array data \eqn{p \times p \times n} in  \eqn{\mathcal{C}_p^+}.
#' @export
#' @importFrom MASS mvrnorm
#'
cor.simulate <- function(n,alpha,sigma,beta = NULL,attC= NULL,p = NULL){
  if(is.null(p)){
    if(!is.null(attC)) {
      p = dim(attC)[1]
    } else {
      stop("Please input the dimension of p")
    }
  }
  m   <- 0.5*p*(p-1)  # dimensions
  q = length(alpha)   # lags
  if(is.list(alpha)) {
    Sca = FALSE
  } else {
    Sca = TRUE
  }
  ###
  ### initial points
  ###
  # We generate the initial data points from [CorrM()]
  C <- array(0,c(p,p,n))
  D = list()
  for(i in 1:(q+1)){
    C[,,i] <- CorrM(p)
    D[[i]]  = D.optimal.position(diag(p),C[,,i])
  }
  ### tangent vectors
  newX <- list()
  for(i in 1:q){
    newX[[i]] <- cor.LogMap(C[,,i],C[,,i+1])
  }
  ### coordinate vectors
  newu <- matrix(0, ncol= m, nrow = n-1)
  for(i in 1:q){
    newu[i,] <- cor.coor_vec(newX[[i]],C[,,i],D[[i]])
  }
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
    ### mean of wrapped Gaussian distribution
    mu <- rep(0,length=m)
    ### Optimal positions
    # sum of multiplication of lag vectors and regression regression
    for(l in 1:q){
      Xtmp <- cor.para_tran(C[,,i-l],C[,,i], newX[[i-l]])
      if(Sca){
        mu <- mu + alpha[l]*cor.coor_vec(Xtmp,C[,,i],D[[i]])
      } else{
        mu <- mu + alpha[[l]]*cor.coor_vec(Xtmp,C[,,i],D[[i]])
      }
    }
    # mean reversion
    if(!is.null(beta)){
      Xtmp = cor.LogMap(C[,,i],attC)
      mu = mu + beta*cor.coor_vec(Xtmp,C[,,i],D[[i]])
    }
    # Euclidean vector
    if(Sca){
      newu[i,] <- rnorm(m,mu,sigma)
    } else {
      newu[i,] <- mvrnorm(1,mu,diag(sigma^2))
    }
    ### orthonormal basis in horizontal space
    S = phi(C[,,i],D[[i]])
    ES = hor.ortho_basis(S)
    # project back to tangent space in quotient manifold
    matV <- matrix(0,p,p)
    for(r in 1:m){
      matV <- matV + newu[i,r] * ES[[r]]
    }
    newX[[i]] <- deri_submer(matV,S)
    # Exponential map
    C[,,i+1] = cor.ExpMap(C[,,i],newX[[i]])
    D[[i+1]] = D.optimal.position(diag(p),C[,,i+1])
  }
  return(C)
}
