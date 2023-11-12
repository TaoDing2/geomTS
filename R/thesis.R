#' Data structure in thesis
#'
#' @param Dat List data
#'
#' @return List data for each patients
#' @export
#'
anti_combinelist <- function(Dat){
  seiinfo <- c(2,2,4,13,4,8,4,4,21,16,2,9,7,2,2,5,2,5)
  npat <- length(seiinfo)
  TSPAT <- list()
  for(pat in 1:npat){
    # No. of seizures
    nsei <- seiinfo[pat]
    # choose seizures for Patient pat
    ind <- c(1:nsei)
    Seizures <- Dat[ind]
    names(Seizures) <- paste("Seizure",1:nsei)
    TSPAT[[pat]] <- Seizures
    # drop off stored seizures
    Dat <- Dat[-ind]
  }
  names(TSPAT) <- paste("Patient ID",1:npat)
  return(TSPAT)
}

