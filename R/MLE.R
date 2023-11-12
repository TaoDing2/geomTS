#' MLE for manifold-adapted model with scalar coefficients
#' @description
#' Model inference by MLE for manifold-adapted model with scalar coefficients.
#'
#' @param lagvec Lagged vector results obtained from functions [euc.coordvec_lags], [spd.coordvec_lags] nor [cor.coordvec_lags].
#' @param maxlag Maximum lag with default value of \code{10}.
#' @param Mean Logical value with default value of \code{FALSE}. When it is \code{TRUE}, we will consider the mean reverting term
#' in our manifold-adapted time series model, vice versa.
#'
#' @return Estimated parameters and model parameters, including \describe{
#' \item{theta}{Estimated parameters of autoregressive terms and mean reverting term (if \code{Mean} is \code{TRUE}).}
#' \item{sigma}{Estimated standard deviation of the white noise in the model.}
#' \item{AIC}{AIC value.}
#' \item{BIC}{BIC value.}
#' \item{loglik}{log-likelihood function value.}
#' }
#'
#' @export
#' @seealso [euc.coordvec_lags], [spd.coordvec_lags], [cor.coordvec_lags].
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
    # model selection and log-likelihood value
    res <- modselectSca(v,lagv,theta,sigma)
    AIC = res$AIC
    BIC = res$BIC
    loglik =  res$loglik
    # store the result
    reslags[[q]] <- list(theta = theta, sigma = sigma,AIC = AIC,BIC = BIC, loglik = loglik)
  }
  names(reslags) <- paste("lag",1:maxlag)
  if(Mean){
    class(reslags) = "MVAR"
  } else{
    class(reslags) = "VAR"
  }
  return(reslags)
}

#' Coefficients in the manifold-adapted model with scalar coefficients
#' @description
#' Estimate coefficients in the manifold-adapted model with scalar coefficients,
#' including coefficients in regressive terms and standard deviation of the white noise.
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#'
#' @return Estimated scalar coefficients, including \describe{
#' \item{theta}{Estimated scalar regressive coefficients.}
#' \item{sigma}{Estimated standard deviation of the white noise.}
#' }
#' @export
#'
EstSca <- function(v,lagv){
  # Estimate alpha
  theta <- EstScaTheta(v,lagv)
  # Estimate sigma
  sigma <- EstScaSig(v,lagv,theta)
  output <- list(theta = theta, sigma = sigma)
  return(output)
}

#' Regressive coefficients in the manifold-adapted model with scalar coefficients
#' @description
#' Compute regressive coefficients in the manifold-adapted model with scalar coefficients.
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors that each of it has the same structure with \code{v}.
#'
#' @return A vector of estimated regressive coefficients.
#' @export
#'
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

#' Standard deviation in the manifold-adapted model with scalar coefficients
#' @description
#' Compute standard deviation of white noise in the manifold-adapted model with scalar coefficients.
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param theta A vector of estimated regressive coefficients from the output of [EstScaTheta()].
#'
#' @return Estimated standard deviation of the white noise in the model.
#' @export
#'
#' @seealso [EstScaTheta()].
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

#' Residual matrix after model fitting with scalar coefficient
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param theta A vector of estimated regressive coefficients from the output of [EstScaTheta()].
#'
#' @return Residual matrix which has the same dimension with  \code{v}.
#' @export
#'
#' @seealso [EstScaTheta()].
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

#' Model selection for manifold-adapted model with scalar coefficients
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param theta A vector of estimated regressive coefficients from the output of [EstScaTheta()].
#' @param sigma Estimated standard deviation of the white noise from the output of function [EstScaSig()].
#'
#' @return AIC , BIC, and log-likelihood function values
#' @export
#' @seealso [EstScaTheta()],[EstScaSig()].
modselectSca <- function(v,lagv,theta,sigma){
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
  # AIC and BIC
  K <- q + 1 # number of parameters
  aic <- 2 * K - 2 * loglik
  bic <- K*log(N) - 2*loglik
  return(list(AIC = aic,BIC = bic, loglik = loglik))
}

#' Print model coefficients in manifold-adapted model with scalar coefficients
#'
#' @param model Results from function [MLE.Sca()]. It is the object of \code{VAR} for the model with autoregressive terms only
#' or \code{MVAR} for the model with autoregressive and mean-reverting term.
#' @param fixed.lag Fixed lag. If it is \code{NULL}, the AIC method is used to suggest model lags.
#' @param modelselection Type of model selection criterion, either \code{AIC} or \code{BIC}. The default value is \code{AIC} and
#' if \code{modelselection} is \code{NULL}.
#' @return A vector of estimated parameters.
#' @export
#'
#' @seealso [MLE.Sca()]
print_Sca <- function(model,fixed.lag = NULL,modelselection = c("AIC","BIC")){
  if(is.null(modelselection)) modelselection = "AIC"
  ### find proper lags
  if(is.null(fixed.lag)){
    abic <- c()
    maxlag <- length(model)
    for(q in 1:maxlag){
      if(modelselection == "AIC") {
        abic[q] <- model[[q]]$AIC
      } else {
        abic[q] <- model[[q]]$BIC
      }
    }
    fixed.lag <- which.min(abic)
  }
  estpar <- c(model[[fixed.lag]]$theta, model[[fixed.lag]]$sigma)
  if(inherits(model,"VAR")){
    names(estpar) <- c(paste("alpha",1:fixed.lag,sep = ""),"sigma")
  } else {
    names(estpar) <- c(paste("alpha",1:fixed.lag,sep = ""),"beta","sigma")
  }
  return(estpar)
}






#' MLE for manifold-adapted model with diagonal coefficients
#' @description
#' Model inference by MLE for manifold-adapted model with diagonal coefficients.
#'
#' @param lagvec Lagged vector results obtained from functions [euc.coordvec_lags], [spd.coordvec_lags] nor [cor.coordvec_lags].
#' @param maxlag Maximum lag with default value of \code{10}.
#' @param Mean Logical value with default value of \code{FALSE}. When it is \code{TRUE}, we will consider the mean reverting term
#' in our manifold-adapted time series model, vice versa.
#'
#' @return Estimated parameters and model parameters, including \describe{
#' \item{Theta}{A list of estimated diagonal coefficients of autoregressive terms and mean reverting term (if \code{Mean} is \code{TRUE}).}
#' \item{Sigma}{A vector of estimated diagonals of standard deviation for the white noise}
#' \item{AIC}{AIC values}
#' \item{BIC}{BIC values}
#' \item{loglik}{log-likelihood function value.}
#' }
#'
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
    res <- modselectDiag(v,lagv,Theta,Sigma)
    AIC = res$AIC
    BIC = res$BIC
    loglik = res$loglik
    # store the result
    reslags[[q]] <- list(Theta = Theta,Sigma = Sigma,AIC = AIC, BIC = BIC, loglik = loglik)
  }
  names(reslags) <- paste("lag",1:maxlag)
  if(Mean){
    class(reslags) = "MVAR"
  } else{
    class(reslags) = "VAR"
  }
  return(reslags)
}


#' Diagonal coefficients in the manifold-adapted model with diagonal coefficients
#' @description
#' Estimate diagonal coefficients in the manifold-adapted model with diagonal coefficients,
#' including coefficients in regressive terms and standard deviation of the white noise.
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#'
#' @return Estimated diagonal coefficients, including \describe{
#' \item{Theta}{A list of estimated diagonal regressive coefficients.}
#' \item{Sigma}{A vector of estimated diagonals of standard deviation for the white noise}
#' }
#' @export
#'
EstDiag <- function(v,lagv){
  # Estimate alpha
  Theta <- EstDiagTheta(v,lagv)
  # Estimate sigma
  Sigma <- EstDiagSig(v,lagv,Theta)
  output <- list(Theta = Theta, Sigma = Sigma)
  return(output)
}

#' Regressive coefficients in the manifold-adapted model with diagonal coefficients
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#'
#' @return A list of estimated diagonal regressive coefficients with the length of \eqn{q}.
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


#' Standard deviation in the manifold-adapted model with diagonal coefficients
#' @description
#' Compute standard deviation of white noise in the manifold-adapted model with diagonal coefficients.
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param Theta A list of estimated diagonal regressive coefficients with the length of \eqn{q} from the output of [EstScaTheta()].
#'
#' @return A vector of diagonal stimated diagonal of Sigma
#' @export
#' @seealso [EstScaTheta()]
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

#' Residual matrix after model fitting with diagonal coefficient
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param Theta A list of estimated diagonal regressive coefficients with the length of \eqn{q} from the output of [EstScaTheta()].
#'
#' @return Residual matrix which has the same dimension with  \code{v}.
#' @export
#' @seealso [EstScaTheta()]
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

#' Model selection values for manifold-adapted model with scalar coefficients
#'
#' @param v A \eqn{n \times m} matrix with \eqn{m} dimensions and \eqn{n} observations.
#' @param lagv A list data with the length of \eqn{q} containing lagged vectors  that each of it has the same structure with \code{v}.
#' @param Theta A list of estimated diagonal regressive coefficients with the length of \eqn{q} from the output of [EstScaTheta()].
#' @param Sigma A vector of estimated diagonals of standard deviation for the white noise from the output of [EstDiagSig].
#'
#' @export
#' @return Model selection and log-likelihood function values
#' @seealso  [EstScaTheta()],[EstDiagSig()].
modselectDiag <- function(v,lagv,Theta,Sigma){
  m <- ncol(lagv[[1]])
  N <- nrow(lagv[[1]]) # No. time points
  q <- length(lagv)       # lags
  Sig <- diag(Sigma^2)
  ### Calculate the Residuals matrix
  ResidMat <- residMat.Diag(v,lagv,Theta)
  ### log-likelihood
  tmp<- sum(sapply(1:N, function(i) t(ResidMat[i,])%*%solve(Sig)%*%ResidMat[i,]))
  loglik <- -0.5 * m * N * log(2*pi) - 0.5*N * log(det(Sig)) - 0.5*tmp
  # AIC and BIC
  K <- (q+1)*m # number of parameters
  aic <- 2 * K - 2 * loglik
  bic <- K*log(N) - 2*loglik
  return(list(AIC = aic, BIC = bic, loglik = loglik))
}

