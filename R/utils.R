#' Diagonal matrix
#' @description
#' Preserve diagonals and discard off-diagonal entries of a square matrix \eqn{M}.
#'
#' @param M A square matrix.
#'
#' @return A square diagonal matrix.
#' @export
#'
#' @examples M = matrix(c(1:9),3,3)
#' Diag(M)
Diag <- function(M){
  return(diag(diag(M)))
}

#' Lower triangular matrix
#' @description Check whether a given squarematrix is lower triangular matrix or not.
#'
#' @param M A square matrix.
#'
#' @return Logical value. Return \code{TRUE} if the given square matrix is a lower triangular matrix. Otherwise, output is \code{FALSE}.
#' @export
#'
#' @examples M = matrix(-c(1:9),3,3)
#' M[upper.tri(M,diag = FALSE)] = 0
#' is.lower.tri(M)
is.lower.tri<-  function (M){
  return(all(M[row(M) > col(M)] != 0))
}

#' Inverse of the square root for a non-singular square matrix
#' @description Compute the inverse of the square root of a non-singular square matrix as \eqn{M^{-1/2}}
#' such that \eqn{MM = S^{-1}} and \eqn{S} is non-singular.
#'
#' @param M A non-singular square matrix.
#'
#' @export
#' @import expm
#'
#' @seealso [sqrtm]
#' @references \insertRef{golub2013matrix}{geomTS}
#'
#' @examples M = CovM(5)
#' sol.sqrtm(M)
sol.sqrtm <- function(M){
  return(solve(sqrtm(M)))
}

#' Trace of a square matrix
#' @description Compute the trace of a square matrix \eqn{M}.
#'
#' @param M A square matrix.
#'
#' @return A numerical value. It is the sum of diagonals for a given square matrix \eqn{M}
#' @export
#'
#' @examples tr(matrix(1:9,3,3))
tr <- function(M){
  return(sum(Diag(M)))
}

#' Symmetrization operator
#' @description Compute the Symmetrizatin operator:  \deqn{Sym(M) = \frac{1}{2}(M+M^T) }.
#'
#' @param M A square matrix.
#'
#' @return A symmetric matrix.
#' @export
#'
#' @examples M = matrix(c(1:9),3,3)
#' Sym(M)
Sym <- function(M){
  return(0.5*(M + t(M)))
}

#' Indices of diagonal positions from a square matrix
#' @description Extract indices of diagonal positions from a square matrix. We number each position from a square matrix.
#' For example a \eqn{3 \times 3} matrix, there are indices of entries from 1 to 9. Therefore, the indices of diagonal positions is \eqn{(1,5,9)}.
#'
#' @param p Row (or column) number of the square matrix.
#'
#' @export
#'
#' @examples diag_ind(3)
#'
diag_ind <- function(p){
  M <- matrix(0,p,p)
  m <- 0.5*p*(p+1)
  M[which(lower.tri(M,diag = T))] <- c(1:m)
  M = t(M)
  return(diag(M))
}


#' Wrapped Gaussian distribution
#'
#' @description Generate a symmetric positive definite matrix from a wrapped Gaussian distribution.
#'
#' @param S A SPD matrix
#' @param sig A positive value. It is the standard deviation of normal distribution with zero mean.
#'
#' @return A SPD matrix
#' @export
#' @importFrom  MASS mvrnorm
#'
#' @details The wrapped Gaussian distribution represents a bump of density around the point \eqn{S} with the degree of spread and shape of the density contours determined by the matrix \eqn{\Sigma}.
#' In fact we will work with a related set of distributions obtained by adding a non-zero mean \eqn{w \in \mathbb{R}^m} to the Gaussian distribution. In this case we write
#' \deqn{U \sim N_{T_S\mathcal{M}_p} (w,\Sigma)}
#' to denote the distribution of the random tangent vector \eqn{U}. The resulting distribution of the random point
#' \deqn{X= \mathrm{Exp}_{S}(U)}
#' closely approximates the wrapped Gaussian with zero mean centred at the point \eqn{\mathrm{Exp}_{S} (\sum_i w_iv_i)}
#' when the vector \eqn{w} is sufficiently small, and \eqn{v_i} is the \eqn{i}th orthonormal basis in \eqn{T_{S}\mathcal{M}_p} and \eqn{i = 1,\ldots,m}.
#'
#' @note \eqn{\Sigma = \sigma^2 I_m}.
#' @references \insertRef{mallasto2018wrapped}{geomTS}
#' @seealso [CovM()]
#'
#' @examples wrapGau(diag(5),0.2)
wrapGau <- function(S,sig){
  p = dim(S)[1]
  m = .5*p*(p+1)
  Y = mvrnorm(1,rep(0,m),sig*diag(m))
  E = spd.ortho_basis(S)
  V = matrix(0,p,p)
  for(i in 1:m){
    V = V + Y[i]*E[[i]]
  }
  return(spd.ExpMap(S,V))
}

#' Full-rank covariance matrix
#'
#' @description A easy way to create an arbitrary \eqn{p \times p} covariance matrix with full rank.
#'
#' @param p Row (or column) number of the symmetric matrix.
#'
#' @details \eqn{M \in \mathbb{R}^{p \times p}}  is a \eqn{p \times p} symmetric matrix
#' which elements \eqn{M_{i,j},i,j = 1,...,p} are simulated from a normal distribution \eqn{N(0,0.5^2)},
#' and then full-rank covariance matrix is computed by
#' \deqn{S = \frac{1}{2}(M+M^T)}
#' \deqn{S = SS^T}
#'
#' @return A covariance matrix with full rank. It is a SPD matrix.
#' @export
#'
#' @seealso [wrapGau()]
#' @examples CovM(5)
#'
CovM <- function(p){
  M = matrix(0,p,p)
  while (min(eigen(M)$values) < 10^{-3}) {
    M <- matrix(c(rnorm(p^2,0,0.5)),p,p)
    M <- 0.5*(M+t(M))
    M <- M%*%t(M)
  }
  return(M)
}

#' Full-rank correlation matrix
#' @description A easy way to create an arbitrary \eqn{p \times p} correlation matrix with full rank from a covariance matrix.
#'
#' @param p Row (or column) number of the symmetric matrix.
#'
#' @details It is computed by transforming a covariance matrix from [CovM()] into a correlation matrix.
#'
#' @return A full-rank correlation matrix. It is also a SPD matrix.
#' @export
#'
#' @seealso [CovM()]
#'
#' @examples CorrM(5)
CorrM <- function(p){
  # construct covariance matrix
  S <- CovM(p)
  # convert  covariance matrix into correlation matrix
  C <- cov2cor(S)
  return(C)
}

#' List object to array object
#' @description Transform list object to array object.
#'
#' @param X  List object with length \eqn{n}. Note that each list contains a \eqn{p \times p} square matrix.
#'
#' @return An \eqn{p \times p \times n} array object.
#' @export
#'
#' @seealso [array2list()]
#' @examples X = lapply(1:10, function(i) CovM(5))
#' S = list2array(X)
#'
list2array <- function(X){
  n <- length(X)
  dims <-  dim(X[[1]])
  S <- array(0,c(dims[1],dims[2],n))
  for(i in 1:n){
    S[,,i] <- X[[i]]
  }
  return(S)
}

#' Array object to list object
#' @description Transform array object to list object
#'
#' @param S An \eqn{p \times p \times n} array object.
#'
#' @return List object with length \eqn{n}. Note that each list contains a \eqn{p \times p} square matrix.
#'
#' @export
#'
#' @seealso [list2array()]
#'
#' @examples S = array(0,c(5,5,10))
#' for(i in 1:10){S[,,i] = CovM(5)}
#' X = array2list(S)
#'
#' @examples
#' data(EEG)
#' S = EEG$seizure1
#' Slist = array2list(S)
#' length(Slist)
array2list <- function(S){
  nlist <- dim(S)[3]
  X <- list()
  for(i in 1:nlist){
    X[[i]] <- S[,,i]
  }
  return(X)
}

#' Array object to vector object
#' @description Transform array object into vector object.
#'
#' @param Dat An \eqn{p \times p \times n} array object.
#' @param Diag Logical value with default value of \code{TRUE}, which means the diagonal entries are included. Otherwise, diagonals are excluded.
#'
#' @return A list data contains \eqn{n} vectors. Dimension of vectors is \eqn{\frac{1}{2}p(p-1)} when \code{Diag = FLASE}, and \eqn{\frac{1}{2}p(p+1)} when \code{Diag = TRUE}.
#' The vectors are the lower triangular entries in the  symmetric matrix.
#'
#' @export
#'
#' @examples
#' ### Covariance matrix
#' S = array(0,c(5,5,10))
#' for(i in 1:10){S[,,i] = CovM(5)}
#' array2vec(S,Diag = TRUE)
#' ### Correlation matrix
#' S = array(0,c(5,5,10))
#' for(i in 1:10){S[,,i] = CorrM(5)}
#' array2vec(S,Diag = FALSE)
array2vec <- function(Dat,Diag = TRUE){
  N <- dim(Dat)[3] # time points
  p <- dim(Dat)[1] # order of matrix
  vec <- list()
  for(i in 1:N){
    tmp <- Dat[,,i]
    vec[[i]] <- tmp[lower.tri(tmp,diag = Diag)]
  }
  return(vec)
}

#' List object to matrix object
#' @description Transform list object to matrix object.
#'
#' @param X List object with length \eqn{n}. Note that each list contains a \eqn{p \times p} square matrix.
#' @param byrow Logical values. When it is \code{TRUE}, we append list data by rows, vice versa.
#'
#' @return A stacked matrix with the dimension of \eqn{(n*p) \times p}
#' @export
#'
#' @examples X = lapply(1:10, function(i) CovM(5))
#' mat = list2matrix(X,byrow = TRUE)
#' mat
list2matrix <- function(X,byrow=TRUE){
  n <- length(X)
  M <- X[[1]]
  if(byrow) {
    for(i in 2:n){
      M <- rbind(M,X[[i]])
    }
  } else {
    for(i in 2:n){
      M <- cbind(M,X[[i]])
    }
  }
  return(M)
}

#' Unified list
#' @description Flatten the list data into a single list.
#'
#' @param X List object with the length of \eqn{n}.
#' @export
#'
#' @return A single list data.
#'
#' @examples X1 = lapply(1:5, function(i) CovM(5))
#' X2 = lapply(1:5, function(i) CovM(5))
#' X = list(X1=X1,X2 =X2)
#' Xlist = flattenlist(X)
#' Xlist
flattenlist <- function(X){
  n <- length(X)
  TS <- X[[1]]
  for(i in 2:n){
    TS <- append(TS,X[[i]])
  }
  return(TS)
}

#' Mean of array data
#' @description Compute the mean of array data by the specific dimensions
#'
#' @param M An \eqn{p \times p \times n} array object.
#' @param dim A choice from \code{dim = (1,2,3)} corresponding to array object.
#' For example \code{dim = 3}, it means we compute array mean by the dimension of \code{3}.
#'
#' @return Matrix mean of dimension \code{dim}
#' @export
#'
#' @examples S = array(0,c(5,5,10))
#' for(i in 1:10){S[,,i] = CovM(5)}
#' mS1 = MeanArray(S,dim = 1)
#' mS2 = MeanArray(S,dim = 2)
#' mS3 = MeanArray(S,dim = 3)
MeanArray <- function(M, dim = NULL){
  d1 = dim(M)[1]
  d2 = dim(M)[2]
  d3 = dim(M)[3]
  if(dim == 1){
    sumM = matrix(0,d2,d3)
    for(i in 1:d1){
      sumM = sumM + M[i,,]
    }
    return(sumM/d1)
  } else if (dim == 2){
    sumM = matrix(0,d1,d3)
    for(i in 1:d2){
      sumM = sumM + M[,i,]
    }
    return(sumM/d2)
  } else {
    sumM = matrix(0,d1,d2)
    for(i in 1:d3){
      sumM = sumM + M[,,i]
    }
    return(sumM/d3)
  }
}
