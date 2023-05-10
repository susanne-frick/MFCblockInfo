#' Calculate D-optimality
#'
#' @param infos block information as returned from calc.info.block()
#' @param summed logical, should the information be summed across blocks? defaults to TRUE
#' @param prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#'
#' @return numeric vector of D-optimality for each person
#' @export
#'
#' @examples
calc.d.optimality <- function(infos, prior = NULL, summed = TRUE) {
  # sum across blocks if there are several
  if(isTRUE(summed) & (length(dim(infos[[1]])) > 2)) {
    infos <- lapply(infos, function(ip) matrix(colSums(ip, dims=1), dim(ip)[2], dim(ip)[3]))
  }
  
  if(isFALSE(is.null(prior))) {
    # invert 
    iprior <- solve(prior)
    # expand the prior if there are several blocks (or not summed)
    if(length(dim(infos[[1]])) > 2) iprior <- array(rep(iprior, each = dim(infos[[1]])[1]), dim = dim(infos[[1]]))
    # add the prior
    infos <- lapply(infos, function(ip, p) ip + p, p = iprior)
  }
  #calculate determinant
  if(length(dim(infos[[1]])) > 2) {
    do.call(rbind, lapply(infos, function(ip) apply(ip, 1, det)))
  } else {
    do.call(c, lapply(infos, function(ip) det(ip)))
  }
  
}
