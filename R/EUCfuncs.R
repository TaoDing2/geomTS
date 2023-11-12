################################################################################
###
### Euclidean geometry on symmetric matrix space
###
################################################################################

#' Euclidean distance in symmetric matrix space
#' @description Calculate Euclidean distance for any two data points in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#'
#' @param S1 A symmetric matrix.
#' @param S2 A symmetric matrix.
#'
#' @details
#' Given any two points \eqn{S_1} and \eqn{S_2}  in symmetric matrix space \eqn{\mathrm{Sym}(p)}, we can define the Euclidean (or Frobenius) distance as
#' \deqn{d(S_1,S_2) =  \| S_1 -S_2 \|_F}
#' where \eqn{\|\cdot \|}  is the Frobenius distance.
#'
#' @return A non-negative number.
#' @export
#'
#' @examples p = 5
#' S1 = CovM(p)
#' S2 = CovM(p)
#' euc.metric(S1,S2)
#' euc.metric(S2,S1)
euc.metric <- function(S1,S2){
  V = S1 - S2
  return(norm(V,"F"))
}

#' Euclidean metric in symmetric matrix space
#' @description Compute Euclidean metric for any two tangent vector matrices in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#'
#' @param V1 A tangent vector matrix.
#' @param V2 A tangent vector matrix.
#'
#' @details A Riemannian metric is defined at \eqn{S \in \mathrm{Sym}(p)} for any two tangent vectors \eqn{V,W \in T_S\mathrm{Sym}(p)= \mathrm{Sym}(p)}
#'  as follows: \deqn{ g_S(V,W) = \mathrm{tr}(V^TW)}
#' This implies the inner product is in fact independent of the base point \eqn{S}.
#'
#' @return A numerical distance.
#' @export
#'
#' @examples p = 5
#' S1 = CovM(p)
#' S2 = CovM(p)
#' S3 = CovM(p)
#' # Compute tangent vectors
#' V1 = S2-S1
#' V2 = S3-S1
#' euc.inner_product(V1,V2)
#' euc.inner_product(V2,V1)
euc.inner_product <- function(V1,V2){
  return(tr(V1%*%V2))
}

#' Logarithm map in symmetric matrix space
#' @description Logarithm map in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#' It aims to obtain a tangent vector matrix by computing the velocity from point \eqn{S_1} to point \eqn{S_2} both locating in \eqn{\mathrm{Sym}(p)}.
#'
#' @param S1 A symmetric matrix.
#' @param S2 A symmetric matrix.
#'
#' @return A tangent vector matrix in tangent space \eqn{T_{S_1}\mathrm{Sym}(p)}.
#' @export
#'
#' @note The tangent vector \eqn{V} also lies in symmetric matrix space \eqn{\mathrm{Sym}(p)}. It is the inverse map of exponential map by [euc.ExpMap()].
#' @seealso [euc.ExpMap()]
#' @examples p = 5
#' S1 = CovM(p)
#' S2 = CovM(p)
#' euc.LogMap(S1,S2)
euc.LogMap <- function(S1,S2){
  return(S2-S1)
}

#' Exponential map in symmetric matrix space
#' @description Exponential map in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#' It projects tangent vector matrix \eqn{V \in T_{S}\mathrm{Sym}(p)} at point \eqn{S \in \mathrm{Sym}(p)} back to \eqn{\mathrm{Sym}(p)}.
#'
#' @param S A symmetric matrix.
#' @param V A tangent vector matrix.
#'
#' @return A symmetric matrix. It is the mapped point of \eqn{V} at base point \eqn{S} in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#' @export
#'
#' @note It is the inverse map of logarithm map by [euc.LogMap()].
#' @seealso [euc.LogMap()]
#' @examples p = 5
#' S1 = CovM(p)
#' S2 = CovM(p)
#' # Compute logarithm map at base point \eqn{S1} lying in tangent space \eqn{T_S\mathrm{Sym}(p)}
#' V = euc.LogMap(S1,S2)
#' euc.ExpMap(S1,V)
euc.ExpMap <- function(S,V){
  return(S + V)
}

#' Geodesic in symmetric matrix space
#' @description Compute geodesic \eqn{\gamma(t), t \in [0,1]} in symmetric matrix space  \eqn{\mathrm{Sym}(p)}.
#'
#' @param t  Time. It should be in the interval of \eqn{t \in [0,1]}.
#' @param S1 A symmetric matrix. It is the starting point of the geodesic \eqn{\gamma(t=0)} in  \eqn{\mathrm{Sym}(p)}.
#' @param S2 A symmetric matrix. It is the endpoint of the geodesic \eqn{\gamma(t=1)} in  \eqn{\mathrm{Sym}(p)}.
#'
#' @return A symmetric matrix. It lies on the geodesic \eqn{\gamma(t)}  at time \eqn{t, t \in [0,1]}.
#' @export
#'
#' @details It is defined by \eqn{\mathrm{Exp}_{S_1}(t V)}, where \eqn{V = \mathrm{Log}_{S_1}(S_2)} and \eqn{t \in [0,1]}.
#'
#' @examples p = 5
#' S1 = CovM(p)
#' S2 = CovM(p)
#' euc.geodesic(t = 0.5,S1,S2)
euc.geodesic <- function(t,S1,S2){
  V = euc.LogMap(S1,S2)
  geodesic = S1 + t*V
  return(geodesic)
}

#' Orthonormal basis at the base point of identity matrix in symmetric matrix space
#' @description Construct orthonormal basis in the tangent space \eqn{T_{I_p}\mathrm{Sym}(p)}.
#'
#' @param p Row (or column) number of the symmetric matrix.
#'
#' @return Orthonormal basis in the tangent space \eqn{T_{I_p}\mathrm{Sym}(p)} in a list object with the dimension of \eqn{m= p(p+1)/2}.
#' @export
#'
#' @details
#' Let \eqn{e_i} denote the \eqn{i}th standard basis vector in \eqn{\mathbb{R}^p},
#' where \eqn{e_i = (0,...,1,...,0)_p^T}. For any \eqn{1 \leq i \leq j \leq p},
#' we can construct an orthonormal basis in \eqn{T_{I_p}\mathrm{Sym}(p)} as
#' \deqn{ E_{ij} = e_ie_j^T \; \mathrm{when} \; i = j,}
#' and \deqn{ E_{ij} = \frac{\sqrt{2}}{2}\left(e_ie_j^T + e_je_i^T\right), \; \mathrm{when} \; i \neq j.}
#' To simplify the notation, let \eqn{E_{ij}} denote as \eqn{E_r} where \eqn{r = 1,...,m= p(p+1)/2}.
#' This represents the \eqn{r}th basis vector in the tangent space \eqn{T_{I_p}\mathrm{Sym}(p)}  at the point \eqn{ I_p \in \mathrm{Sym}(p)}.
#'
#' @note According to the property of Euclidean geometry \eqn{\mathrm{Sym}(p)},
#' the orthonormal basis in tangent spaces does not depend on the base point ,
#' i.e., for given any \eqn{S \in \mathrm{Sym}(p)},
#' the orthonormal basis in the tangent space \eqn{T_{S}\mathrm{Sym}(p)} is exactly same as basis in the tangent space \eqn{T_{I_p}\mathrm{Sym}(p)} .
#'
#' @examples euc.ortho_basis(p = 5)
euc.ortho_basis <- function(p){
  E <- list()
  m <- 0.5*p*(p+1)
  num <- 1
  for(i in 1:p){
    for(j in i:p){
      e1 <- rep(0,length= p)
      e2 <- rep(0,length =p)
      e1[i] <- e2[j] <- 1
      if( i == j)
        E[[num]] <- 0.5* (e1%*%t(e2) + e2%*%t(e1))
      else E[[num]] <- 0.5*sqrt(2)*(e1%*%t(e2) + e2%*%t(e1))
      num <- num + 1
    }
  }
  return(E)
}

#' Coordinate vector by vectorizing projection in symmetric matrix space
#' @description Vectorize a tangent vector matrix \eqn{V \in T_S\mathrm{Sym}(p)} into a coordinate vector \eqn{\boldsymbol{v} \in \mathbb{R}^m}. Thereby,
#' the dimension of symmetric matrix space is \eqn{m = p(p+1)/2}. If the symmetric matrix space consists of correlation matrices,
#' the dimension of symmetric matrix space will be \eqn{m = p(p-1)/2}.
#'
#'
#' @param V A tangent vector matrix.
#' @param corMat Logical values with default value \code{FALSE}.
#' If it is \code{TRUE}, symmetric matrix space is consisting of correlation matrix and the dimension of coordinate vector \eqn{\boldsymbol{v}} is
#'  \eqn{m =p(p-1)/2}. Otherwise, \eqn{m = p(p+1)/2}.
#'
#' @return A coordinate vector in \eqn{\mathbb{R}^m}
#' @export
#'
#' @seealso [euc.Euc2Sym()]
#' @examples S1 = CovM(5)
#' S2 = CovM(5)
#' V = euc.LogMap(S1,S2)
#' euc.coor_vec(V)
euc.coor_vec <- function(V,corMat = NULL){
  if( is.null(corMat)) corMat <- FALSE
  p = ncol(V)
  m <- 0.5* p * (p + 1)
  E <- euc.ortho_basis(p) # orthonormal basis on T_{I_p}
  v <- sapply(1:m, function(i) euc.inner_product(V,E[[i]]))
  if(corMat){
    return(v[-diag_ind(p)])
  } else {
    return(v)
  }
}

#' Mapping back the coordinate vector to the data point in symmetric matrix space
#' @description Map back the coordinate vector to the tangent space \eqn{T_{I_p}\mathrm{Sym}(p)}
#' and then do exponential map at base point \eqn{S \in \mathrm{Sym}(p)}.
#'
#' @param v A coordinate vector in \eqn{\mathbb{R}^m}.
#' @param S A symmetric matrix in \eqn{\mathrm{Sym}(p)}.
#' @param corMat Logical values with default value \code{FALSE}. If it is \code{TRUE},
#' symmetric matrix space is consisting of correlation matrix and the dimension of coordinate vector \eqn{\boldsymbol{v}} is
#'  \eqn{m =p(p-1)/2}. Otherwise, \eqn{m = p(p+1)/2}.
#'
#' @return A symmetric matrix in \eqn{\mathrm{Sym}(p)}
#' @export
#'
#' @seealso [euc.coor_vec()]
#' @examples S1 = CovM(5)
#' S2 = CovM(5)
#' v = euc.coor_vec(S2-S1)
#' euc.Euc2Sym(v,S1)
#'
#' # Correlation matrix
#' C1 = CorrM(5);C2 = CorrM(5)
#' u = euc.coor_vec(C2-C1,corMat = TRUE)
#' euc.Euc2Sym(u,C1,corMat = TRUE)
euc.Euc2Sym <- function(v,S,corMat = NULL){
  if( is.null(corMat)) corMat <- FALSE
  # dimension
  m = length(v)
  p = ncol(S)
  # basis in tangent space
  E = euc.ortho_basis(p)
  if(corMat) {
    E = E[-diag_ind(p)]
  }
  ### Tangent vector
  V = matrix(0,p,p)
  for(i in 1:m){
    V = V + v[i]*E[[i]]
  }
  ### point in manifold
  return(euc.ExpMap(S,V))
}


#' Euclidean mean in symmetric matrix space
#' @description Compute Euclidean mean  in symmetric matrix space \eqn{\mathrm{Sym}(p)}
#'
#' @param S Array object with the dimension of \eqn{p \times p \times n}. It includes \eqn{n} observations of \eqn{p \times p} symmetric matrices.
#'
#' @return A symmetric matrix in \eqn{\mathrm{Sym}(p)}. It is the centre point of dataset \eqn{S}.
#' @export
#'
#' @examples S = lapply(1:15, function(i) CovM(5))
#' euc.mean(list2array(S))
euc.mean <- function(S){
  return(MeanArray(S,3))
}
