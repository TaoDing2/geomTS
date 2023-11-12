#' MLE for VAR model with scalar coefficients
#'
#' @param lagvec lagged results obtained from functions [spd.coordvec_lags] or [cor.coordvec_lags]
#' @param maxlag maximum lags with default 10
#' @param Mean Logical values. True for involving Mean reversion, vice versa.
#'
#' @return estimated parameters and AIC values for different lags
#' @export
#'
# examples S = lapply(1:21, function(i) CovM(5))
# S = list2array(S)
# FSM = spd.mean_Sturm(S, MaxIt = 100)
# lagvec = spd.coordvec_lags(S,FSM,maxlag = 3)
# mod = MLE.Sca(lagvec,maxlag = 3)
MLE.Sca <- function(lagvec,maxlag = NULL,Mean = NULL){
  if(is.null(maxlag)) maxlag = 10
  if(is.null(Mean)) Mean = FALSE
  ### Estimate parameters for different lags
  reslags <- list()
  for(q in 1:maxlag){
    # vectors
    v = lagvec$v[-(1:q),]
    # Obtain the lagged vectors
    lagv <- lagvec$lagv[[q]]
    if(Mean){
      lagv[[q+1]] = lagvec$mv[-c(1:q),]
    }
    # obtain estimated parameters
    estpar <- EstSca(v,lagv)
    theta <- estpar$theta
    sigma <- estpar$sigma
    # AIC and log-likelihood value
    AICres<- AICSca(v,lagv,theta,sigma)
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

#' Estimate Parameters in the model with scalar coefficients
#'
#' @param v vectors
#' @param lagv lagged vectors
#'
#' @return estimated scalar coefficients
#' @export
#'
#' @seealso The whole function could find in function [MLE.Sca()]
EstSca <- function(v,lagv){
  # Estimate alpha
  theta <- EstScaTheta(v,lagv)
  # Estimate sigma
  sigma <- EstScaSig(v,lagv,theta)
  output <- list(theta = theta, sigma = sigma)
  return(output)
}

#' Estimate alpha
#'
#' @param v vectors
#' @param lagv lagged vectors
#'
#' @return estimated alpha
#' @export
#'
#' @seealso The whole function could find in function [MLE.Sca()] and [EstSca()]
EstScaTheta <- function(v,lagv){
  N <- nrow(lagv[[1]])
  m <- ncol(lagv[[1]])
  q <- length(lagv)
  if(q ==1){                # Z matrix with m by q
    Y <- sum(sapply(1:N, function(i) v[i,]%*%lagv[[1]][i,]))
    Z <- sum(sapply(1:N, function(i) lagv[[1]][i,]%*%lagv[[1]][i,]))
    theta <- Y/Z
  } else {
    Y <- c()
    Z <- matrix(0,q,q)
    for(l in 1:q){
      Y[l] <- sum(sapply(1:N, function(i) v[i,]%*%lagv[[l]][i,]))
      for(o in 1:q){
        Z[l,o] <- sum(sapply(1:N, function(i) lagv[[l]][i,]%*%lagv[[o]][i,]))
      }
    }
    theta <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
  }
  return(theta)
}

#' estimate Sigma in scalar model
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param theta estimated alpha
#'
#' @return estimated sigma
#' @export
#'
#' @seealso The whole function could find in function [MLE.Sca()],[EstSca()],[EstScaTheta()]
EstScaSig <- function(v,lagv,theta){
  N <- nrow(lagv[[1]])
  m <- ncol(lagv[[1]])
  # Calculate the Residuals matrix
  ResidMat <- residMat.Sca(v,lagv,theta)
  # Calculate residual sum of squares (RSS)
  RSS <- 0
  for(i in 1:N){
    RSS <- RSS + ResidMat[i,]%*%ResidMat[i,]
  }
  sig <- sqrt(RSS/(m*N))
  return(sig)
}

#' Compute residual matrix for scalar coefficient model
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param theta estimated alpha
#'
#' @return residual matrix
#' @export
#'
residMat.Sca <- function(v,lagv,theta){
  N <- nrow(lagv[[1]])
  m <- ncol(lagv[[1]])
  q <- length(lagv)
  # Calculate the Residuals matrix
  ResidMat <- matrix(0,nrow= N,ncol= m)
  for(i in 1:N){
    lsum <- rep(0,length = m)
    for(l in 1:q){
      lsum <- lsum + theta[l]*lagv[[l]][i,]
    }
    ResidMat[i,] <- v[i,] - lsum
  }
  return(ResidMat)
}

#' Model selection by AIC
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param theta estimated alpha
#' @param sigma estimated sigma
#'
#' @return AIC values and log-likelihood
#' @export
#'
AICSca <- function(v,lagv,theta,sigma){
  N <- nrow(lagv[[1]])
  m <- ncol(lagv[[1]])
  q <- length(lagv)
  # Calculate the Residuals matrix
  ResidMat <- residMat.Sca(v,lagv,theta)
  # Calculate residual sum of squares (RSS)
  RSS <- 0
  for(i in 1:N){
    RSS <- RSS + ResidMat[i,]%*%ResidMat[i,]
  }
  # log-likelihood value
  loglik <- - 0.5 * m * N * log(2*pi) - m * N * log(sigma) - 0.5 * (1/sigma^2) * RSS
  # AIC and corrected AIC
  K <- q+1 # number of parameters
  aic <- 2 * K - 2 * loglik
  return(list(AIC = aic, loglik = loglik))
}


#' Outputs of model with scalar coefficients
#'
#' @param model results from function [MLE.Sca()]. The object could be either "MVAR" or "VAR"
#' @param fixed.lag choose the results by a fixed lag
#'
#' @return estimated parameters
#' @export
#'
print_Sca <- function(model,fixed.lag = NULL){
  ### find proper lags
  if(is.null(fixed.lag)){
    aic <- c()
    maxlag <- length(model)
    for(q in 1:maxlag){
      aic[q] <- model[[q]]$AIC
    }
    fixed.lag <- which.min(aic)
  }
  estpar <- c(model[[fixed.lag]]$theta, model[[fixed.lag]]$sigma)
  if(inherits(model,"VAR")){
    names(estpar) <- c(paste("alpha",1:fixed.lag,sep = ""),"sigma")
  } else {
    names(estpar) <- c(paste("alpha",1:fixed.lag,sep = ""),"beta","sigma")
  }
  return(estpar)
}

################################################################################
####
#### MLE for the model with diagonal matrix coefficients
####
#' MLE for the model with diagonal matrix coefficients
#'
#' @param lagvec lagged vectors
#' @param maxlag maximum lags
#' @param Mean Logical values. True for involving Mean reversion, vice versa.
#'
#' @return  estimated parameters and AIC values for different lags
#' @export
#'
MLE.Diag <- function(lagvec,maxlag = NULL, Mean = NULL){
  if(is.null(maxlag)) maxlag = 10
  if(is.null(Mean)) Mean = FALSE
  ### Estimate parameters for different lags
  reslags <- list()
  for(q in 1:maxlag){
    # vectors
    v = lagvec$v[-(1:q),]
    # Obtain the lagged vectors
    lagv <- lagvec$lagv[[q]]
    if(Mean){
      lagv[[q+1]] = lagvec$mv[-c(1:q),]
    }
    # obtain estimated parameters
    estpar <- EstDiag(v,lagv)
    Theta <- estpar$Theta
    Sigma <- estpar$Sigma
    # AIC
    AICres <- AICDiag(v,lagv,Theta,Sigma)
    AIC = AICres$AIC
    loglik = AICres$loglik
    # store the result
    reslags[[q]] <- list(Theta = Theta,Sigma = Sigma,AIC = AIC, loglik = loglik)
  }
  names(reslags) <- paste("lag",1:maxlag)
  if(Mean){
    class(reslags) = "MVAR"
  } else{
    class(reslags) = "VAR"
  }
  return(reslags)
}


#' Estimate Parameters in diagonal matrix VAR model
#'
#' @param v vectors
#' @param lagv lagged vectors
#'
#' @export
EstDiag <- function(v,lagv){
  # Estimate alpha
  Theta <- EstDiagTheta(v,lagv)
  # Estimate sigma
  Sigma <- EstDiagSig(v,lagv,Theta)
  output <- list(Theta = Theta, Sigma = Sigma)
  return(output)
}

#' Estimate diagonal matrices A in VAR model
#'
#' @param v vectors
#' @param lagv lagged vectors
#'
#' @return estimated diagonals of coefficients
#' @export
#'
EstDiagTheta <- function(v,lagv){
  m <- ncol(lagv[[1]])
  N <- nrow(lagv[[1]]) # No. time points
  q <- length(lagv)       # lags
  ###
  ### construct matrix of Y
  ###
  Y <- c()                          # store Y
  #  calculate u_i * v_i
  u <- list()
  for(l in 1:q){
    tmp <- rep(0,length=m)
    for(i in 1:N){
      tmp <- tmp + diag(lagv[[l]][i,])%*%v[i,]
    }
    u[[l]] <- tmp
  }
  if(q==1){
    Y <- u[[1]]
  } else {
    Y <- u[[1]]
    for(l in 2:q){
      Y <- rbind(Y,u[[l]])
    }
  }
  ###
  ### construct matrix of Z
  ###
  Dia <- list()
  num <- 1
  for(l in 1:q){
    for (o in 1:q) {
      tmp <- matrix(0,m,m)
      for(i in 1:N){
        tmp <- tmp + diag(lagv[[l]][i,])%*%diag(lagv[[o]][i,])
      }
      Dia[[num]] <- tmp
      num <- num +1
    }
  }
  num <- 1
  Z <- matrix(0,ncol= q*m, nrow= q*m)
  for(l in 1:q){
    for(o in 1:q){
      Z[((l-1)*m + 1):(l*m),((o-1)*m + 1):(o*m)] <- Dia[[num]]
      num <- num +1
    }
  }
  ###
  ### estimate A*
  ###
  alp <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
  Theta <- list()
  for(l in 1:q){
    Theta[[l]] <- alp[((l-1)*m + 1):(l*m)]
  }
  return(Theta)
}


#' Estimate diagonals of Sigma in VAR model
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param Theta estimated diagonals of coefficients in VAR models
#'
#' @return estimated diagonal of Sigma

#' @export
EstDiagSig <- function(v,lagv,Theta){
  m <- ncol(lagv[[1]])
  N <- nrow(lagv[[1]]) # No. time points
  q <- length(lagv)       # lags
  ###
  ### Estimate Sigma
  ###
  ResiM <- matrix(0,ncol = m,nrow=N)
  for(i in 1:N){
    sumM <- rep(0,length=m)
    for(l in 1:q){
      sumM <- sumM + diag(lagv[[l]][i,])%*%Theta[[l]]
    }
    ResiM[i,] <- (v[i,] - sumM)^2
  }
  Sig <- c()
  tmp <- colSums(ResiM)/N
  for(i in 1:m){
    Sig[i] <- sqrt(tmp[i])
  }
  return(Sig)
}

#' Compute residual matrix for diagonal matrix coefficient
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param Theta estimated alpha
#'
#' @return residual matrix
#' @export
#'
residMat.Diag <- function(v,lagv,Theta){
  m <- ncol(lagv[[1]])
  N <- nrow(lagv[[1]]) # No. time points
  q <- length(lagv)       # lags
  ###
  ### Estimate Sigma
  ###
  ResidMat <- matrix(0,ncol = m,nrow=N)
  for(i in 1:N){
    sumM <- rep(0,length=m)
    for(l in 1:q){
      sumM <- sumM + diag(lagv[[l]][i,])%*%Theta[[l]]
    }
    ResidMat[i,] <- v[i,] - sumM
  }
  return(ResidMat)
}

#' AIC and log-likelihood
#'
#' @param v vectors
#' @param lagv lagged vectors
#' @param Theta estimated A
#' @param Sigma estimated Sigma
#'
#' @export
AICDiag <- function(v,lagv,Theta,Sigma){
  m <- ncol(lagv[[1]])
  N <- nrow(lagv[[1]]) # No. time points
  q <- length(lagv)       # lags
  Sig <- diag(Sigma^2)
  ### Calculate the Residuals matrix
  ResidMat <- residMat.Diag(v,lagv,Theta)
  ### log-likelihood
  tmp<- sum(sapply(1:N, function(i) t(ResidMat[i,])%*%solve(Sig)%*%ResidMat[i,]))
  loglik <- -0.5 * m * N * log(2*pi) - 0.5*N * log(det(Sig)) - 0.5*tmp
  # AIC and corrected AIC
  K <- (q+1)*m # number of parameters
  aic <- 2 * K - 2 * loglik
  return(list(AIC = aic, loglik = loglik))
}

