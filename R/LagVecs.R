#' Coordinate vectors with lags in symmetric matrix space
#' @description Compute coordinate vectors in  \eqn{\mathbb{R}^m} including observed tangent vectors,
#' mean reversion vectors and lagged vectors in symmetric matrix space  \eqn{\mathrm{Sym}(p)}.
#'
#' @param S  An \eqn{p \times p \times n} array data consisting of symmetric matrices.
#' @param maxlag Maximum lag. The default value is 10.
#' @param FSM A symmetric matrix. It could be the Frechet sample mean \eqn{S} or an attract point.
#' If this value is missing, the mean reversion coordinate vector will not be considered in the model.
#' @param corMat Logical values. If it is \code{TRUE}, the dimension of coordinate vector is
#' \eqn{p(p+1)/2 } , otherwise, \eqn{p(p-1)/2}. The default value is \code{FALSE}.
#'
#' @export
#'
#' @return A set of coordinate vectors, including lagged vectors, observed vectors,
#'  mean reversion vectors (If FSM is not \code{NULL}), and original matrices.
#'
euc.coordvec_lags <- function(S,FSM = NULL,maxlag = NULL,corMat = NULL){
  # default  maximum lag = 10
  if( is.null(maxlag)) maxlag <- 10
  if( is.null(corMat)) corMat <- FALSE
  N <- dim(S)[3]-1 # No. of time points
  p <- dim(S)[1]
  if(corMat){
    m <- 0.5* p * (p - 1)
  } else {
    m <- 0.5* p * (p + 1)
  }
  ### stop when the maximum lags are equal or less than total of observations
  if(N <= maxlag){
    stop("Input a smaller maxlag")
  }
  # Euclidean vectors
  v = matrix(0,N,m)
  for(i in 1:N){
    v[i,] = euc.coor_vec(euc.LogMap(S[,,i],S[,,i+1]),corMat)
  }
  ###
  ### mean vectors
  ###
  if(!missing(FSM)) {
    mv <- matrix(0,nrow = N, ncol = m)
    for(i in 1:N){
      mv[i,] = euc.coor_vec(euc.LogMap(S[,,i],FSM),corMat)
    }
  }
  ###
  ### lagged coordinate vectors
  ###
  lagvs = list()
  for(q in 1:maxlag){
    ### lagged vectors
    lagv = list()
    for(i in 1:q){
      lagv[[i]] = v[(q-i+1):(N-i),]
    }
    lagvs[[q]] = lagv
  }

  ###
  ### return values
  ###
  if(missing(FSM)) {
    lagvec = list(lagv = lagvs, v = v,M = S)
  } else {
    lagvec = list(lagv = lagvs,mv = mv, v=v,M = S)
  }
  return(lagvec)
}


#' @title  Coordinate vectors with lags in affine invariant geometry
#'
#' @description Compute coordinate vectors in  \eqn{\mathbb{R}^{p(p+1)/2}} including observed tangent vectors,
#' mean reversion vectors and lagged vectors in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#'
#' @param S  An \eqn{p \times p \times n} array data consisting of SPD matrices.
#' @param FSM A SPD matrix. It could be the Frechet sample mean \eqn{S} or an attract point.
#' If this value is missing, the mean reversion coordinate vector will not be considered in the model.
#' @param maxlag Maximum lags. The default value is 10.
#'
#' @export
#'
#' @return A set of coordinate vectors, including lagged vectors, observed vectors,
#'  mean reversion vectors (If FSM is not \code{NULL}), and original matrices.
#'
#' @seealso [spd.meanvecs()],[spd.paraV()],[spd.lagvecs()]
spd.coordvec_lags <- function(S,FSM =NULL,maxlag = NULL){
  # default  maximum lag = 10
  if( is.null(maxlag)) maxlag <- 10
  N <- dim(S)[3] - 1 # No. of time points
  p <- dim(S)[1]
  m <- 0.5*p*(p+1)
  ### stop when the maximum lags are equal or less than total of observations
  if(N <= maxlag){
    stop("Input a smaller maxlag")
  }
  ###
  ### Tangent vectors and coordinate vectors
  ###
  cat("Computing tangent vectors and coordinate vectors","\n")
  V <- list()
  v <- matrix(0,ncol= m,nrow=N)
  for(i in 1:N){
    # Tangent vectors
    V[[i]] <- spd.LogMap(S[,,i],S[,,i+1])
    # coordinate vectors
    v[i,] <- spd.coor_vec(V[[i]],S[,,i])
  }
  ###
  ### Mean reversion vectros and coordinate vectors
  ###
  if(!missing(FSM)){
    cat("Computeing mean reversion vectors","\n")
    mv = spd.meanvecs(S,FSM)
  }
  ###
  ### parallel transported tangent vectors
  ###
  cat("Computeing parallel transport vectors","\n")
  paraV <- spd.paraV(maxlag,S,V)
  ###
  ### lagged coordinate vectors
  ###
  cat("Computing lagged coordinate vectors","\n")
  lagvs <- list()
  cat("lags = ")
  for(q in 1:maxlag){
    if(q == maxlag){
      cat(q,"\n")
    } else {
      cat(q,",",sep = "")
    }
    # Calculate the lagged vectors
    lagvs[[q]] <- spd.lagvecs(q,S,paraV)
  }
  ###
  ### return values
  ###
  if(missing(FSM)) {
    lagvec = list(lagv = lagvs, v = v,M = S)
  } else {
    lagvec = list(lagv = lagvs,mv = mv, v = v,M = S)
  }
  return(lagvec)
}

#' Coordinate vectors of mean reversion in affine invariant geometry
#' @description
#' Compute coordinate vectors at the centre point or attract point in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#'
#' @details
#' Let \eqn{\tilde{S}} be Frechet sample mean or attract point. We can compute the tangent vectors for all data points into tangent space \eqn{T_{\tilde{S}}\mathcal{S}_p^+} by
#' \deqn{V_i = \mathrm{Log}_{\tilde{S}}(S_i), i = 1,\ldots,n.}
#' Subsequently, we have the coordinate vectors by vectorization \eqn{Vec(V_i)} by [spd.coor_vec()].
#'
#' @param S An \eqn{p \times p \times n} array data consisting of SPD matrices.
#' @param FSM A SPD matrix. It could be the Frechet sample mean of data \eqn{S} or a given attract point.
#'
#' @return A matrix containing coordinate vectors of mean reversion with the row of observations and the column of the dimension \eqn{m}.
#' @export
#'
#' @seealso [spd.coordvec_lags()]
#'
spd.meanvecs <- function(S,FSM){
  p = dim(S)[1]
  m = 0.5*p*(p+1)
  N <- dim(S)[3]-1
  mv <- matrix(0,ncol = m,nrow = N)
  for(i in 1:N){
    Vtmp <- spd.LogMap(S[,,i],FSM)
    mv[i,] <- spd.coor_vec(Vtmp,S[,,i])
  }
  return(mv)
}

#' Parallel transport vectors with maximum lags in affine invariant geometry
#' @description Compute lagged tagent vector matrix by parallel transport in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#'
#'
#' @param maxlag A numerical value which represents the maximum lags involved
#' @param S An (\eqn{p \times p \times n}) array data consisting of SPD matrices.
#' @param V A (\eqn{n-1}) length list data consisting of tangent vectors
#'
#' @details Given any original data point \eqn{S} in \eqn{\mathcal{S}_p^+} and corresponding
#' tangent vectors \eqn{V}, we compute all parallel transport vectors with the maximum lags \eqn{maxlag}.
#' For example, \code{maxlag = 3}, the parallel transported vector could be written as
#' \eqn{V_{4,1} = \mathcal{P}_{S_1,S_4}(V_1)}, \eqn{V_{4,2} = \mathcal{P}_{S_2,S_4}(V_2)}, and \eqn{V_{4,3} = \mathcal{P}_{S_3,S_4}(V_3)}.
#' The basic structure (e.g maxlag = 3) stored in outputs is
#'
#' \eqn{\mathcal{P}_{S_1,S_2}(V_1)};
#'
#' \eqn{\mathcal{P}_{S_2,S_3}(V_2)}, \eqn{\mathcal{P}_{S_1,S_3}(V_1)};
#'
#' \eqn{\mathcal{P}_{S_3,S_4}(V_3)}, \eqn{\mathcal{P}_{S_2,S_4}(V_2)} ,\eqn{\mathcal{P}_{S_1,S_4}(V_1)};
#'
#' \eqn{\mathcal{P}_{S_4,S_5}(V_4)}, \eqn{\mathcal{P}_{S_3,S_5}(V_3)} ,\eqn{\mathcal{P}_{S_2,S_5}(V_2)}.
#'
#' @return A list object for the parallel transported coordinate vectors
#' @export
#' @seealso [spd.lagvecs()]
spd.paraV <- function(maxlag,S,V){
  N = length(V)
  paraV <- list()
  for(i in 2:N){
    tmp_paraV = list()
    # parallel transport maximum lag is employed
    if(i > maxlag){
      for(j in (i-1):(i-maxlag)){
        tmp_paraV[[i-j]] = spd.para_tran(S[,,j],S[,,i],V[[j]])
      }
    } else{
      for(j in (i-1):1){
        tmp_paraV[[i-j]] = spd.para_tran(S[,,j],S[,,i],V[[j]])
      }
    }
    paraV[[i-1]] = tmp_paraV
  }
  return(paraV)
}

#' Lagged coordinate vectors from parallel transport in affine invariant geometry
#' @description Compute lagged coordinate vectors from parallel transport and arrange the structure of coordinate vectors.
#'
#' @param q A numerical value of lags.
#' @param S An (\eqn{p \times p \times n}) array data consisting of SPD matrices.
#' @param paraV Parallel transported tangent vectors from the output of [spd.paraV()].
#'
#' @return Lagged \eqn{q} coordinate vectors with \eqn{q} lists.
#' Each list is a matrix with the row of translated coordinate vectors and the column of dimensions of coordinate vectors.
#'
#' @export
#'
#' @seealso [spd.coordvec_lags()],[spd.paraV()]
#'
spd.lagvecs <- function(q,S,paraV){
  N = length(paraV)
  p = dim(S)[1]
  m = 0.5*p*(p+1)
  lagv = list() # list coordinate vectors with lag q
  for(l in 1:q){
    vtmp <- matrix(0,ncol= m,nrow= (N-q+1))
    for(i in q:N){
      tmpV = paraV[[i]][[l]]
      vtmp[i-q+1,] = spd.coor_vec(tmpV,S[,,i+1])
    }
    lagv[[l]] = vtmp
  }
  return(lagv)
}


#' Coordinate vectors with lags in quotient geometry
#' @description
#' Compute coordinate vectors in  \eqn{\mathbb{R}^{p(p-1)/2}} including observed tangent vectors,
#' mean reversion vectors and lagged vectors in quotient geometry \eqn{\mathcal{C}_p^+}.
#'
#' @param C An (\eqn{p \times p \times n}) array data consisting of correlation matrices in \eqn{\mathcal{C}_p^+}.
#' @param FSM A correlation matrix. It could be the Frechet sample mean S or an attract point.
#' If this value is missing, the mean reversion coordinate vector will not be considered in the model.
#' @param maxlag Maximum lags. The default value is 10.
#'
#' @return A set of coordinate vectors, including lagged vectors, observed vectors,
#'  mean reversion vectors (If FSM is not \code{NULL}), and original matrices.
#' @export
#' @seealso [cor.meanvecs()],[cor.paraX()],[cor.lagvecs()]
cor.coordvec_lags <- function(C,FSM = NULL,maxlag = NULL){
  # default  maximum lag = 10
  if( is.null(maxlag)) maxlag <- 10
  N <- dim(C)[3] -1 # No. of time points
  p <- dim(C)[1] # dimensions of correlation matrix
  I <- diag(p) # identity matrix
  m <- 0.5* p * (p - 1)       # dimensions of coordinate vectors
  ### stop when the maximum lags are equal or less than total of observations
  if(N <= maxlag){
    stop("Input a smaller maxlag")
  }
  ###
  ### Finding optimal positions with respect to \eqn{I_p}
  ###
  cat("Finding optimal positions","\n")
  D = list()
  for(i in 1:N){
    D[[i]] = D.optimal.position(I,C[,,i])
  }
  ###
  ### Tangent vectors and coordinate vectors
  ###
  cat("Computing tangent vectors","\n")
  X <- list()
  u <- matrix(0,nrow= N,ncol= m)
  for(i in 1:N){
    # Tangent vectors
    X[[i]] <- cor.LogMap(C[,,i],C[,,i+1])
    # coordinate vectors
    u[i,] <- cor.coor_vec(X[[i]],C[,,i],D[[i]])
  }
  ###
  ### Mean reversion vectors and coordinate vectors
  ###
  if(!missing(FSM)){
    cat("Computeing mean reversion vectors","\n")
    mu = cor.meanvecs(C,FSM,D)
  }
  ###
  ### parallel transported tangent vectors
  ###
  cat("Computeing parallel transport vectors","\n")
  paraX <- cor.paraX(maxlag,C,X)
  ###
  ### lagged coordinate vectors
  ###
  cat("Computing lagged coordinate vectors","\n")
  lagvs <- list()
  cat("lags = ")
  for(q in 1:maxlag){
    if(q == maxlag){
      cat(q,"\n")
    } else {
      cat(q,",",sep ="")
    }
    # Calculate the lagged vectors
    lagvs[[q]] <- cor.lagvecs(q,C,paraX,D)
  }
  ###
  ### return values
  ###
  if(missing(FSM)) {
    lagvec = list(lagv = lagvs, v = u,M = C)
  } else {
    lagvec = list(lagv = lagvs,mv = mu, v = u,M = C)
  }
  return(lagvec)
}

#' Coordinate vectors of mean reversion in quotient geometry
#' @description
#' Compute coordinate vectors at the centre point or attract point in quotient geometry  \eqn{\mathcal{C}_p^+}.
#'
#' @param C An (\eqn{p \times p \times n}) array data consisting of correlation matrices in \eqn{\mathcal{C}_p^+}.
#' @param FSM A correlation matrix. It could be Frechet sample mean or  a given attract point.
#' @param D Diagonal matrices in \eqn{\mathcal{D}_p^+}.
#' It ensures \eqn{DC_iD} and \eqn{I_p} are in optimal positions in \eqn{\mathcal{S}_p^+}, \eqn{i = 1,\ldots,n}.
#'
#' @return A matrix containing coordinate vectors of mean reversion with the row of observations and the column of the dimension \eqn{m}
#' @export
#'
#' @seealso [cor.coordvec_lags()]
cor.meanvecs <- function(C,FSM,D){
  p = dim(C)[1]
  m = 0.5*p*(p-1)
  N <- dim(C)[3]-1
  mv <- matrix(0,ncol = m,nrow=N)
  for(i in 1:N){
    Xtmp = cor.LogMap(C[,,i],FSM)
    mv[i,] = cor.coor_vec(Xtmp,C[,,i],D[[i]])
  }
  return(mv)
}

#' Parallel transport vectors with maximum lags in quotient geometry
#' @description Compute lagged tagent vector matrix by parallel transport in affine invariant geometry \eqn{\mathcal{C}_p^+} .
#'
#' @param maxlag A numerical value which represents the maximum lags involved.
#' @param C An (\eqn{p \times p \times n}) array data consisting of correlation matrices in \eqn{\mathcal{C}_p^+}.
#' @param X A \eqn{n-1} length list data consisting of tangent vectors.
#'
#' @details Given any original data point \eqn{C} in \eqn{\mathcal{C}_p^+} and corresponding
#' tangent vectors \eqn{X}, we compute all parallel transport vectors with the maximum lags \eqn{maxlag}.
#' For example, \code{maxlag = 3}, the parallel transported vector could be written as
#' \eqn{X_{4,1} = \mathcal{P}_{C_1,C_4}(X_1)}, \eqn{X_{4,2} = \mathcal{P}_{C_2,C_4}(X_2)}, and \eqn{X_{4,3} = \mathcal{P}_{C_3,S_4}(C_3)}.
#' The basic structure (e.g maxlag = 3) stored in outputs is
#' \eqn{\mathcal{P}_{C_1,C_2}(X_1)};
#'
#' \eqn{\mathcal{P}_{C_2,C_3}(X_2)}, \eqn{\mathcal{P}_{C_1,C_3}(X_1)};
#'
#' \eqn{\mathcal{P}_{C_3,C_4}(X_3)}, \eqn{\mathcal{P}_{C_2,C_4}(X_2)} ,\eqn{\mathcal{P}_{C_1,C_4}(X_1)};
#'
#' \eqn{\mathcal{P}_{C_4,C_5}(X_4)}, \eqn{\mathcal{P}_{C_3,C_5}(X_3)} ,\eqn{\mathcal{P}_{C_2,C_5}(X_2)}
#'
#' @return The parallel transported vectors.
#' @export
#'
#' @seealso [cor.lagvecs()]
cor.paraX <- function(maxlag,C,X){
  N = length(X)
  paraX <- list()
  for(i in 2:N){
    tmp_paraX = list()
    # parallel transport maximum lag is employed
    if(i > maxlag){
      for(j in (i-1):(i-maxlag)){
        tmp_paraX[[i-j]] = cor.para_tran(C[,,j],C[,,i],X[[j]])
      }
    } else{
      for(j in (i-1):1){
        tmp_paraX[[i-j]] = cor.para_tran(C[,,j],C[,,i],X[[j]])
      }
    }
    paraX[[i-1]] = tmp_paraX
  }
  return(paraX)
}


#' Lagged coordinate vectors from parallel transport in quotient geometry
#' @description
#' Compute lagged coordinate vectors from parallel transport and arrange the structure of coordinate vectors in quotient geometry \eqn{\mathcal{C}_p^+}
#'
#' @param q A numerical value of lags.
#' @param C An (\eqn{p \times p \times n}) array data consisting of correlation matrices in \eqn{\mathcal{C}_p^+}.
#' @param paraX  Parallel transported tangent vectors from the putput of [cor.paraX()].
#' @param D Diagonal matrices in \eqn{\mathcal{D}_p^+}.
#' It ensures \eqn{DC_iD} and \eqn{I_p} are in optimal positions in \eqn{\mathcal{S}_p^+}, \eqn{i = 1,\ldots,n}.
#'
#' @return Lagged \eqn{q} coordinate vectors with \eqn{q} lists.
#' Each list is a matrix with the row of translated coordinate vectors and the column of dimensions of coordinate vectors.
#' @export
#'
#' @seealso [cor.coordvec_lags()],[cor.paraX()]
cor.lagvecs <- function(q,C,paraX,D){
  N <- length(paraX)
  p = dim(C)[1]
  m = 0.5*p*(p-1)
  lagv <- list() # list coordinate vectors with lag q
  for(l in 1:q){
    utmp <- matrix(0,ncol=m,nrow=(N-q+1))
    for(i in q:N){
      tmpX = paraX[[i]][[l]]
      utmp[i-q+1,] = cor.coor_vec(tmpX,C[,,i+1],D[[i+1]])
    }
    lagv[[l]] = utmp
  }
  return(lagv)
}





