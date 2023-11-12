#' OLS for model estimation
#'
#' @param v observed coordinated vectors
#' @param mv Mean vectors
#' @param lagv Lagged vectors
#' @param intercept Logical values. If it is TRUE, the model consider the drift
#'
#' @return the regression results which same as [lm()]
#' @export
#' @import utils
#' @import stats
#'
OLS.Sca <- function(v, mv = NULL, lagv = NULL, intercept = NULL){
  ### check the default values
  if(is.null(intercept)) intercept = FALSE
  if(!is.null(lagv)) {
    AR = TRUE
    q = length(lagv)
  } else AR = FALSE
  if(!is.null(mv)) MR = TRUE
  else MR = FALSE
  ### observations
  N = dim(v)[1]
  ### variables
  m = dim(v)[2]
  ### consider AR and MR
  if(MR){
    ### check should we consider Mean reverting regression model
    if(AR){
      ### X
      X1 =lapply(1:q, function(l) stack(as.data.frame(t(lagv[[l]])))$values)
      X2 = list(stack(as.data.frame(t(mv)))$values)
      X = append(X1,X2)
      X = as.data.frame(X)
      names(X) = c(paste("L",1:q,sep = ""),"M")
    } else {
      ### X
      X = list(stack(as.data.frame(t(mv)))$values)
      X = as.data.frame(X)
      names(X) = c("M")
    }
  } else {
    if(AR) {
      ### X
      X =lapply(1:q, function(l) stack(as.data.frame(t(lagv[[l]])))$values)
      X = as.data.frame(X)
      names(X) = paste("L",1:q,sep = "")
    } else {
      print("You should input one of mv and lagv at least")
    }
  }
  if(intercept) {
    indX = as.data.frame(kronecker(rep(1,N),diag(m)))
    X = append(indX,X)
  }
  ### Y
  Y = stack(as.data.frame(t(v)))$values
  Y = as.data.frame(Y)
  ### Combine X and Y
  dat = append(Y,X)
  ### Linear regression
  reg =  lm(Y~ -1 + .,data = dat)
  return(reg)
}
