#' convert information to SEs
#'
#' across persons
#'
#' @param infos list of arrays with block information as returned from calc.info.block()
#' @param summed logical, should the information be summed across blocks? defaults to TRUE
#' @param var.out logical, return variances instead of SEs? defaults to FALSE
#' @param prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#'
#' @return if summed across bocks: matrix: rows=persons, columns=traits,
#' otherwise: list with entries for each person: matrix of SEs: rows = blocks, columns = traits
#' @export
#'
info2se <- function(infos, summed=TRUE, var.out=FALSE, prior=NULL) {
  
  #SEs from infos (within a person)
  .info2se <- function(infos, summed=TRUE, var.out, prior) {
    if (isTRUE(summed)) infos <- matrix(colSums(infos, dims=1), dim(infos)[2], dim(infos)[3])
    
    #several blocks
    if (length(dim(infos))==3) {
      if(is.null(prior)) {
        ses <- t(apply(infos, 1, function(info) diag(MASS::ginv(info))))
      } else {
        ses <- t(apply(infos, 1, function(info) diag(MASS::ginv(info + solve(prior)))))
      }
      if(isFALSE(var.out)) {
        ses <- t(apply(ses, 1, function(se) sqrt(se)))
      }
      
      #summed / only one block
    } else {
      if(is.null(prior)) {
      ses <- diag(MASS::ginv(infos))
    } else {
      ses <- diag(MASS::ginv(infos + solve(prior)))
    }
      if(isFALSE(var.out)) {
        ses <- sqrt(ses)
      }
    }
    return(ses)
  }
  
  if (is.list(infos)) {
    ses <- lapply(infos, .info2se, summed=summed, var.out=var.out, prior=prior)
  } else {
    ses <- .info2se(infos, summed=summed, var.out=var.out, prior=prior)
  }
  
  if(isTRUE(summed) | length(ses) == 1) {
    return(do.call(rbind, ses))
  } else {
    sesre <- matrix(do.call(c, ses), nrow = length(ses), byrow = TRUE)
    return(array(sesre, dim = c(length(ses), dim(ses[[1]])[1], dim(ses[[1]])[2])))
  }
}
