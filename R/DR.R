#' Dimensional reduction via maximizing variation
#' @description
#' Dimensional reduction for full rank covariance (or correlation) matrices via
#' PCA method.
#'
#' @param X An \eqn{p \times p \times n} array data.
#' @param Y An \eqn{p \times p \times n} array data with the same structure with \eqn{X}. If it is \code{NULL}, we will do PCA method for \eqn{X}, otherwise,
#' we will reduce dimensions of \eqn{Y} by multiplying the eigenvector matrix from \eqn{X}.
#' @param N0 Reduced dimension.
#' @param CORmat Logical values. Default value is \code{FALSE},
#' the return data will be covariance matrices, otherwise, correlation matrices.
#'
#' @return An \eqn{N_0 \times N_0 \times n} array data.
#' @export
#'
#' @details
#' Computing the full \eqn{q \times q} sample covariance matrices \eqn{\{S'_i, i = 1,\ldots,n\}},
#'  we then seek linear combinations of channels \eqn{u \in \mathbb{R}^q} which maximize \eqn{\frac{1}{n}\sum_{i}u^TS_i'u},
#'  where \eqn{(\cdot)^T} denotes the transpose.
#'  In other words, we find the eigenvectors \eqn{u_1,\ldots,u_p} of \eqn{\frac{1}{n}\sum_{i}S_i'} with the largest\eqn{p} eigenvalues.
#'  The time series of covariance (or correlation) matrices \eqn{S_i} (or \eqn{C_i}) defined as the sample covariances (correlation) of these linear combinations:
#'  \deqn{ S_i  =  \mathrm{Cov}{\{Uz_{(i-1)f + 1},\ldots,Uz_{if}\}} = US_i'U^T }
#'  or
#'  \deqn{C_i  = \mathrm{Cor}{ \{Uz_{(i-1)f + 1},\ldots,Uz_{if}\} }}
#'  where the \eqn{p \times q} matrix \eqn{U} has rows \eqn{u_1^T,\ldots,u_p^T}.
#'
#'
#' @note
#' \eqn{X} and \eqn{Y} should be the covariance matrices.
#' When we reduce the dimension of correlation matrices,
#' we tend to do PCA on covariance matrices and then convert it to correlation matrices.
#' @examples
#' data(EEG)
#' # Covariance matrix
#' eegcov = EEG$seizure1
#' resu = DR_PCA(eegcov,N0 =15,CORmat =FALSE)
#' # Correlation matrix using common eigenvectors
#' eegcov2 = EEG$interictal1
#' resu2 = DR_PCA(X = eegcov,N0 = 15,Y = eegcov2,CORmat = TRUE)
DR_PCA <- function(X,N0,Y = NULL,CORmat = FALSE){
  # No of observations
  N = dim(X)[3]
  # means in Euclidean space
  MX = MeanArray(X,dim = 3)
  ###
  ###  PCA for dimensionality reduction
  ###
  # ### Eigenvalues
  # D = eigen(MX)$values
  ### Eigenvectors
  EV = eigen(MX)$vectors
  # dim-reduction
  W = EV[,1:N0]
  ###
  ### Reduced data
  ###
  # If Y is Null, it means we do PCA to X independently, otherwise, Y will use
  # eigenvector matrix obtained from X
  if(is.null(Y)) {
    if(CORmat){
      Mat = lapply(1:N, function(j) cov2cor(t(W)%*%X[,,j]%*%W))
    } else {
      Mat = lapply(1:N, function(j) t(W)%*%X[,,j]%*%W)
    }
  } else {
    if(CORmat){
      Mat = lapply(1:N, function(j) cov2cor(t(W)%*%Y[,,j]%*%W))
    } else {
      Mat = lapply(1:N, function(j) t(W)%*%Y[,,j]%*%W)
    }
  }
  return(list2array(Mat))
}



#' Dimensionality reduction via minimizing redundancy
#' @description
#' Dimensional reduction for full rank covariance (or correlation) matrices by
#' choosing the greatest mean of minimum eigenvalues (GME) via greedy algorithm.
#'
#' @param X An \eqn{p \times p \times n} array data.
#'
#' @return The descending order with choosing pairs and corresponding GMEs
#' @export
#'
#'
#' @details
#' Suppose we have a set of channels \eqn{\mathfrak{C}}, and we let \eqn{S_i'(\mathfrak{C})}
#' (or \eqn{C_i'(\mathfrak{C})}) denote the restriction of the full covariance (or correlation)
#' matrix \eqn{S_i'} (or \eqn{C'_i}) to \eqn{\mathfrak{C}}. Also define
#'   \deqn{\omega(\mathfrak{C})  = \frac{1}{n}\sum_{i}\min \sigma[S_i'(\mathfrak{C})] }
#'   or \deqn{\omega(\mathfrak{C})  = \frac{1}{n}\sum_{i}\min \sigma[C_i'(\mathfrak{C})]}
#'  where \eqn{\sigma(A)} denotes the set of eigenvalues of \eqn{A \in  \mathrm{Sym}(p)}.
#'  A greedy algorithm was used to construct sets \eqn{\mathfrak{C}} which maximized \eqn{\omega(\mathfrak{C})} for a fixed value of \eqn{p}.
#'  More specifically, we first consider every set \eqn{\mathfrak{C}} consisting of a pair of channels to find the optimal set \eqn{\mathfrak{C}_2}.
#'   We then consider every set \eqn{\mathfrak{C}} containing \eqn{\mathfrak{C}_2}
#'   and one additional channel, to find the optimal set \eqn{\mathfrak{C}_3}, and so on.
#'
#' @examples
#' data(EEG)
#' eegcov = EEG$seizure1
#' DR_GME(eegcov)
DR_GME <- function(X){
  p <- dim(X)[1] # No. of channels
  N <- dim(X)[3] # No. of time points
  g.m.ev <- c() # store the greatest mean minimum eigenvalues
  ###
  ###  Step 1: Find two "best" channels
  ###
  ### combinations of pairs of channels
  inipairs <- combn(p,2)
  npairs <- dim(inipairs)[2]
  mevals <- rep(0,npairs)
  for(i in 1:npairs){
    channels <- inipairs[,i]
    lambda = c()
    for(t in 1:N){
      S0 = X[channels,channels,t]
      lambda[t] = min(eigen(S0)$values)
    }
    # store means of minimum eigenvalues
    mevals[i] <- mean(lambda)
  }
  ### Find the greatest mean minimum eigenvalues
  g.m.ev[1]<- max(mevals)
  # The pair of channels which has the greatest mean of minimum eigenvalues
  pair <- inipairs[,which(mevals == max(mevals),arr.ind = TRUE)]
  ###
  ### Step 2: Add the channels based on greedy algorithm
  ###
  for(itr in 2:(p-1)){
    #    cat("Finding ",itr+1,"th channel","\n")
    ind <- c(1:p)
    ind <- ind[-pair]
    mevals <- rep(0,length=length(ind))
    for(i in 1:length(ind)){
      channels <- c(pair, ind[i]) # extract matrices with chal X chal diemnsion
      lambda <- c()
      for(t in 1:N){
        S0 = X[channels,channels,t]
        lambda[t] = min(eigen(S0)$values)
      }
      # store means of minimum eigenvalues
      mevals[i] <- mean(lambda)
    }
    # The greatest mean minimum eigenvalues with increasing channels
    g.m.ev[itr] <- max(mevals)
    # The combination of channels which has the greatest mean minimum eigenvalues
    pair <- c(pair, ind[which(mevals == max(mevals))]) # increasing channels added
  }
  Pairs <- list(pairs = pair,g.m.ev = g.m.ev)
  return(Pairs)
}

#' Dimensionality reduction by choosing the highest minimum eigenvalues via greedy algorithm
#'
#' @param X An \eqn{p \times p \times n} array data.
#'
#' @return The descending order of HME with choosing pairs and corresponding HMEs
#' @export
#' @examples
#' data(EEG)
#' eegcov = EEG$seizure1
#' DR_HME(eegcov)
DR_HME <- function(X){
  p <- dim(X)[1] # No. of channels
  N <- dim(X)[3] # No. of time points
  h.m.ev <- c() # store the highest minimum eigenvalues
  ###
  ###  Step 1: Find two "best" channels
  ###
  ### combinations of pair of channels
  inipairs <- combn(p,2)
  npairs <- dim(inipairs)[2]
  minevals <- rep(0,npairs)
  for(i in 1:npairs){
    channels <- inipairs[,i]
    tmp <- min(eigen(X[channels,channels,1])$values)
    for(t in 2:N){
      S = X[channels,channels,t]
      lambda = min(eigen(S)$values)
      if(lambda < tmp){
        tmp = lambda
      }
    }
    # store minimum eigenvalues
    minevals[i] <- tmp
  }
  ### Find the highest minimum eigenvalues
  h.m.ev[1]<- max(minevals)
  # The pair of channels which has the highest minimum eigenvalues
  pair <- inipairs[,which(minevals == max(minevals),arr.ind = TRUE)]
  ###
  ### Step 2: Add the channels based on greedy algorithm
  ###
  for(itr in 2:(p-1)){
    #    cat("Finding ",itr+1,"th channel","\n")
    ind <- c(1:p)
    ind <- ind[-pair]
    minevals <- rep(0,length=length(ind))
    for(i in 1:length(ind)){
      channels <- c(pair, ind[i]) # extract matrices with chal X chal diemnsion
      tmp <- min(eigen( X[channels,channels,1])$values)
      for(t in 2:N){
        S = X[channels,channels,t]
        lambda = min(eigen(S)$values)
        if(lambda < tmp){
          tmp = lambda
        }
      }
      # store minimum eigenvalues
      minevals[i] <- tmp
    }
    # The highest minimum eigenvalues with increasing channels
    h.m.ev[itr] <- max(minevals)
    # The combination of channels which has the highest minimum eigenvalues
    pair <- c(pair, ind[which(minevals == max(minevals))]) # increasing channels added

  }
  Pairs <- list(pairs = pair,h.m.ev = h.m.ev)
  return(Pairs)
}
