
#' Weighted MLE
#'
#' @param LAGVEC Lagged coordinate vectors from the outputs of [euc.coordvec_lags()], [spd.coordvec_lags()], and [cor.coordvec_lags()]
#' @param maxlag A numerical value of maximum lag
#' @param Mean Logical value with default \code{FALSE} that the regression model does not consider mean reversion term.
#'
#' @return A "MVAR" or "VAR" object containing \describe{
#' \item{theta}{Estimated regression coefficients}
#' \item{sigma}{Estimated s.d.}
#' \item{AIC}{AIC values}
#' \item{liglik}{log-likelihood values}
#' }
#'
#' @export
#'
wMLE.Sca <- function(LAGVEC,maxlag = NULL, Mean = NULL){
  if(is.null(maxlag)) maxlag = 5
  if(is.null(Mean)) Mean = FALSE
  # No of datasets
  nsei <- length(LAGVEC)
  ### weights
  omega <- sapply(1:nsei, function(sei) nrow(LAGVEC[[sei]]$v))
  omega <- omega/sum(omega)
  ### Estimate parameters for different lags
  reslags <- list()
  for(q in 1:maxlag){
    # obtain estimated parameters
    theta <- EstScaTheta4all(q,omega,LAGVEC,Mean)
    sigma = EstScaSig4all(q,omega,LAGVEC,theta, Mean)
    # AIC and log-likelihood value
    AICres <- AICSca4for(q,omega,LAGVEC,theta,sigma, Mean)
    AIC = AICres$AIC
    loglik =  AICres$loglik
    # store the result
    reslags[[q]] <- list(theta = theta, sigma = sigma,AIC = AIC,liglik = loglik)
  }
  names(reslags) <- paste("lag",1:maxlag)
  if(Mean){
    class(reslags) = "MVAR"
  } else{
    class(reslags) = "VAR"
  }
  return(reslags)
}


#' Estimated regression coefficients
#'
#' @param q A numerical value with lags
#' @param omega A vector of weights
#' @param LAGVEC Lagged coordinate vectors from the outputs of [euc.coordvec_lags()], [spd.coordvec_lags()], and [cor.coordvec_lags()]
#' @param Mean Logical value with default \code{FALSE} that the regression model does not consider mean reversion term.
#'
#' @return Estimated regression coefficients
#' @export
EstScaTheta4all <- function(q,omega,LAGVEC,Mean = NULL){
  if(is.null(Mean)) Mean = FALSE
  # number of seizures
  nsei = length(LAGVEC)
  Ysei = list()
  Zsei =list()
  ### compute the Y and Z matrices for each seizure
  for(sei in 1:nsei){
    v = LAGVEC[[sei]]$v[-c(1:q),]
    lagv = LAGVEC[[sei]]$lagv[[q]]
    if(Mean){
      mv = LAGVEC[[sei]]$mv[-c(1:q),]
      lagv[[q+1]] = mv
      nq = q + 1
    } else {
      nq = q
    }
    N = nrow(v)
    if(q ==1){                # Z matrix with m by q
      Y <- sum(sapply(1:N, function(i) v[i,]%*%lagv[[1]][i,]))
      Z <- sum(sapply(1:N, function(i) lagv[[1]][i,]%*%lagv[[1]][i,]))
    } else {
      Y <- c()
      Z <- matrix(0,nq,nq)
      for(l in 1:nq){
        Y[l] <- sum(sapply(1:N, function(i) v[i,]%*%lagv[[l]][i,]))
        for(o in 1:nq){
          Z[l,o] <- sum(sapply(1:N, function(i) lagv[[l]][i,]%*%lagv[[o]][i,]))
        }
      }
    }
    Ysei[[sei]] = omega[sei]* Y
    Zsei[[sei]] = omega[sei]* Z
  }
  ### Weighted  Y and Z matrices
  Y = Reduce(`+`, Ysei)
  Z = Reduce(`+`, Zsei)
  if(q==1){
    theta = Y/Z
  } else {
    theta = solve(t(Z)%*%Z)%*%t(Z)%*%Y
  }
  return(theta)
}


#' Estimated s.d.
#'
#' @param q A numerical value with lags
#' @param omega A vector of weights
#' @param LAGVEC Lagged coordinate vectors from the outputs of [euc.coordvec_lags()], [spd.coordvec_lags()], and [cor.coordvec_lags()]
#' @param theta Estimated regression coefficients
#' @param Mean Logical value with default \code{FALSE} that the regression model does not consider mean reversion term.
#'
#' @return Estimated s.d.
#' @export
EstScaSig4all <- function(q,omega,LAGVEC,theta, Mean = NULL){
  if(is.null(Mean)) Mean = FALSE
  # number of seizures
  nsei = length(LAGVEC)
  RSSsei = 0
  Coef = 0
  for(sei in 1:nsei){
    v = LAGVEC[[sei]]$v[-c(1:q),]
    N = nrow(v)
    m = ncol(v)
    lagv = LAGVEC[[sei]]$lagv[[q]]
    if(Mean){
      mv = LAGVEC[[sei]]$mv[-c(1:q),]
      lagv[[q+1]] = mv
    }
    # Calculate the Residuals matrix
    ResidMat <- residMat.Sca(v,lagv,theta)
    # Calculate residual sum of squares (RSS)
    RSS <- 0
    for(i in 1:N){
      RSS <- RSS + ResidMat[i,]%*%ResidMat[i,]
    }
    RSSsei  = RSSsei +  omega[sei]*RSS
    Coef = Coef + omega[sei]* m*N
  }
  return(sqrt(RSSsei/Coef))
}


#' AIC and log-likelihood values
#'
#' @param q A numerical value with lags
#' @param omega A vector of weights
#' @param LAGVEC Lagged coordinate vectors from the outputs of [euc.coordvec_lags()], [spd.coordvec_lags()], and [cor.coordvec_lags()]
#' @param theta Estimated regression coefficients
#' @param sigma Estimated s.d.
#' @param Mean Logical value with default \code{FALSE} that the regression model does not consider mean reversion term.
#'
#' @return AIC and log-likelihood values
#' @export
#'
AICSca4for <- function(q,omega,LAGVEC,theta,sigma, Mean = NULL){
  if(is.null(Mean)) Mean = FALSE
  # number of seizures
  nsei = length(LAGVEC)
  loglik = c()
  for(sei in 1:nsei){
    v = LAGVEC[[sei]]$v[-c(1:q),]
    N = nrow(v)
    m = ncol(v)
    lagv = LAGVEC[[sei]]$lagv[[q]]
    if(Mean){
      mv = LAGVEC[[sei]]$mv[-c(1:q),]
      lagv[[q+1]] = mv
    }
    # Calculate the Residuals matrix
    ResidMat <- residMat.Sca(v,lagv,theta)
    # Calculate residual sum of squares (RSS)
    RSS <- 0
    for(i in 1:N){
      RSS <- RSS + ResidMat[i,]%*%ResidMat[i,]
    }
    # log-likelihood value
    loglik[sei] <- - 0.5 * m * N * log(2*pi) - m * N * log(sigma) - 0.5 * (1/sigma^2) * RSS
  }
  # AIC and corrected AIC
  K <- q  + 1 # number of parameters
  aic <- 2 * K - 2 * loglik*omega
  return(list(AIC = sum(aic), loglik = sum(loglik*omega)))
}
