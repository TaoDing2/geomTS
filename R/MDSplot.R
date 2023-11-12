#' Distance matrix using Euclidean metric
#' @description
#' Compute distance matrix of symmetric matrices in symmetric matrix space \eqn{\mathrm{Sym}(p)}.
#'
#' @param S An \eqn{p \times p \times n} array data consisting of symmetric matrices.
#'
#' @return Distance object as [dist()].
#' @export
#'
#' @seealso [euc.dist()]
#' @examples
#' data(EEG)
#' redEEG = DR_PCA(EEG$seizure1,N0 = 15)
#' euc.dist(redEEG)
euc.dist <- function(S){
  N = dim(S)[3]
  d <- matrix(0,N,N)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      d[j,i] = euc.metric(S[,,i],S[,,j])
    }
  }
  return(as.dist(d))
}

#' Distance matrix using affine invariant metric
#' @description
#' Compute distance matrix of SPD matrices in affine invariant geometry \eqn{\mathcal{S}_p^+}.
#'
#' @param S An \eqn{p \times p \times n} array data consisting of SPD matrices.
#' @param parallel Logical value with default \code{FALSE}.When it is \code{TRUE},
#' we use parallel computation to compute distance matrix.
#'
#' @return Distance object as [dist()].
#' @export
#'
#' @import doParallel
#' @import foreach
#'
#' @seealso [spd.dist()]
#' @examples
#' data(EEG)
#' redEEG = DR_PCA(EEG$seizure1,N0 = 15)
#' spd.dist(redEEG)
spd.dist <- function(S,parallel = NULL){
  N = dim(S)[3]
  if(is.null(parallel)) parallel = FALSE
  if(parallel) {
    # combination of matrix indices
    ind = combn(N,2)
    n = ncol(ind)
    # No. of cores used
    registerDoParallel(detectCores()-2)
    dc <- foreach(i = 1:n)%dopar%{
      return(spd.metric(S[,,ind[1,i]],S[,,ind[2,i]]))
    }
    d = matrix(0,N,N)
    d[lower.tri(d,diag = FALSE)] = unlist(dc)
  } else {
    d <- matrix(0,N,N)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        d[j,i] = spd.metric(S[,,i],S[,,j])
      }
    }
  }
  return( as.dist(d))
}

#' Distance matrix using quotient metric
#' @description
#' Compute distance matrix of full rank correlation matrices in quotient geometry \eqn{\mathcal{C}_p^+}.
#'
#' @param S An \eqn{p \times p \times n} array data consisting of full rank correlation matrices.
#' @param parallel Logical value with default \code{FALSE}.When it is \code{TRUE},
#' we use parallel computation to compute distance matrix.
#' @return Distance object as [dist()].
#' @export
#'
#' @import doParallel
#' @import foreach
#'
#' @seealso [cor.dist()]
#' @examples
#' S = lapply(1:10, function(i) CorrM(5))
#' S = list2array(S)
#' cor.dist(S)
cor.dist <- function(S,parallel = NULL){
  N = dim(S)[3]
  if(is.null(parallel)) parallel = FALSE
  if(parallel) {
    # combination of matrix indices
    ind = combn(N,2)
    n = ncol(ind)
    # No. of cores used
    registerDoParallel(detectCores()-2)
    dc <- foreach(i = 1:n)%dopar%{
      return(cor.metric(S[,,ind[1,i]],S[,,ind[2,i]]))
    }
    d = matrix(0,N,N)
    d[lower.tri(d,diag = FALSE)] = unlist(dc)
  } else {
    d <- matrix(0,N,N)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        d[j,i] = cor.metric(S[,,i],S[,,j])
      }
    }
  }
  return( as.dist(d))
}

#' Multidimensional scaling plot
#'
#' @param dist An object of class "dist"
#' @param titlename An title name of multidimensional scaling (MDS) plot. If it is \code{NULL},
#' the title of ggplot will not be shown in the output.
#' @param type MDS methods implemented. The default values is metric MDS \code{type = "mMDS"}.
#' Users can also choose classical MDS \code{type = "cMDS"} and non-metric MDS \code{type = "nMDS"}.
#' @param nc Number of colours. If it is \code{NULL}, we will use 3 gradient colours for scatter plot.
#'
#' @return A ggplot output
#'
#' @export
#'
#' @importFrom scales percent
#' @importFrom colorspace rainbow_hcl
#' @importFrom smacof smacofSym
#' @importFrom MASS isoMDS
#' @import ggplot2
#'
MDSplot <- function(dist, titlename = NULL,type = NULL,nc = NULL){
  if(is.null(type)) type = "mMDS"
  if(type == "cMDS"){
    ### classical MDS scale
    mds = cmdscale(dist,k=2, eig = T)
    x = mds$points[,1]
    y = mds$points[,2]
    N = length(x)
    Time = c(1:N)
    # proportion of variance captured by the fitted vectors
    varc = percent(mds$GOF[2],accuracy = 0.01)
    df = data.frame(x = x,y = y,Time = Time)
  } else if (type == "mMDS") {
    ### metric MDS scale
    mds = smacofSym(dist,ndim =2)
    x = mds$init[,1]
    y = mds$init[,2]
    N = length(x)
    Time = c(1:N)
    # proportion of variance captured by the fitted vectors
    varc = percent(mds$stress,accuracy = 0.01)
    df = data.frame(x = x,y = y,Time = Time)
  } else if (type == "nMDS"){
    ### non-metric MDS scale
    mds = isoMDS(dist,y = cmdscale(dist,k=2),k=2)
    x = mds$points[,1]
    y = mds$points[,2]
    N = length(x)
    Time = c(1:N)
    # proportion of variance captured by the fitted vectors
    varc = percent(mds$stress/100,accuracy = 0.01)
    df = data.frame(x = x,y = y,Time = Time)
  } else {
    print("Please use a correct MDS method!")
  }
  if(is.null(nc)) nc = 3
  # colour bar
  col.rng <- rainbow_hcl(nc,c = 80,l = 50)
  fig = ggplot(df,aes(x,y,colour = Time)) + geom_point(size = 2) +
    geom_path(alpha = 0.3) +
    scale_colour_gradientn(colours = col.rng)
  if(is.null(titlename)) {
    return(fig)
    }
  else {
    return(fig + ggtitle(paste(titlename,paste("(",varc,")",sep = ""),sep=" " )))
  }
}




