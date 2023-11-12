#' \eqn{R^2} for checking the goodness of fit
#' @description Compute \eqn{R^2} for evaluating the goodness of model fitting.
#'
#' @param v Coordinate vectors of observed data.
#' @param lagvs Coordinate vectors of lagged vectors.
#' @param theta Estimated parameters. If it is a list, the model with diagonal matrix coefficients is used.
#' Otherwise, it is the model with scalar coefficients.
#'
#' @return \eqn{R^2}
#' @export
#'
squaredR <- function(v,lagvs,theta){
  q =  length(lagvs)
  m <- dim(v)[2]
  N <- dim(v)[1]
  # Calculate the Residuals matrix
  if(is.list(theta)  ){
    ResidMat <- residMat.Diag(v,lagvs,theta)
  } else{
    ResidMat <- residMat.Sca(v,lagvs,theta)
  }
  # residual sum of squares (RSS)
  RSS <- 0
  for(i in 1:N){
    RSS <- RSS + t(ResidMat[i,])%*%ResidMat[i,]
  }
  # total sum of squares (TSS)
  mv = colMeans(v)
  TSS = 0
  for(i in 1:N){
    TSS <- TSS + t(v[i,] - mv)%*%(v[i,] - mv )
  }
  ### R^2
  return(1- RSS/TSS)
}


#' @title Compute 95% confidence interval w.r.t. VAR model with scalar coefficients
#'
#' @param alpha Estimated alpha
#' @param sigma Estimated sigma
#' @param lagv lagged vectors
#'
#' @export
#' @return Confidence matrix with \eqn{m \times 2} where \eqn{m} is the sum of coefficients of
#'
CalCIs <- function(alpha,sigma,lagv){
  q <- length(lagv)
  N <- nrow(lagv[[1]])
  m <- ncol(lagv[[1]])
  IM <- matrix(0, ncol = q+1,nrow= q+1)
  for(i in 1:q){
    for(j in 1:q){
      IM[i,j] <- as.vector(lagv[[i]])%*%as.vector(lagv[[j]])/(sigma^2) ### i,j
    }
  }
  IM[q+1,q+1] <- 2*m*N/(sigma^2)   ### alpha sigma
  II <- solve(IM)
  CI <- matrix(0,ncol=2,nrow = q+1)
  for(i in 1:q){
    CI[i,] <- c(alpha[i] - 1.96*sqrt(II[i,i]), alpha[i] + 1.96*sqrt(II[i,i]))
  }
  CI[q+1,] <- c(sigma - 1.96*sqrt(II[q+1,q+1]), sigma + 1.96*sqrt(II[q+1,q+1]))
  return(CI)
}


#' Asymptotic covariance of maximum log-likelihood estimator
#' @description Its is computed by the Hessian matrix operator and solved by
#' the inverse of Fisher information matrix.
#'
#' @param sigma Estimated sigma with respect to noise term
#' @param lagv Lagged vectors
#'
#' @return Asymptotic covariance
#' @export
#'
asympCOV <- function(sigma,lagv){
  q <- length(lagv)
  m <- nrow(lagv[[1]])
  n <- ncol(lagv[[1]])
  IM <- matrix(0, ncol = q+1,nrow= q+1)
  for(i in 1:q){
    for(j in 1:q){
      IM[i,j] <- as.vector(lagv[[i]])%*%as.vector(lagv[[j]])/(sigma^2)
    }
  }
  IM[q+1,q+1] <- 2*m*n/(sigma^2)
  return( solve(IM))
}


