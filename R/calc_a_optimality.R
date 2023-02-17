#' calculate A-optimality
#'
#' @param infos block information as returned from calc.block.infos
#' @param summed logical, should the information be summed across blocks? defaults to TRUE
#' @param prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#'
#' @return numeric vector of A-optimality for each person
#' @export
#'
#' @examples
calc.a.optimality <- function(infos, prior=NULL, summed=TRUE) {
  ses <- info2se(infos, var.out=T, prior=prior, summed=summed)
  if(length(dim(ses)) == 3) {
    rowSums(ses, dims = 2)    
  } else if (length(dim(ses)) == 1){
    return(ses) 
  } else {
    rowSums(ses, na.rm=T)
  }
}
