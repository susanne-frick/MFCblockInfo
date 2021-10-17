#' convert information to SEs
#'
#' across persons
#'
#' @param infos list of arrays with block information as returned from calc.info.block()
#' @param summed logical, should the information be summed across blocks? defaults to TRUE
#' @param var.out logical, return variances instead of SEs? defaults to FALSE
#' @param s.prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#'
#' @return if summed across bocks: matrix: rows=persons, columns=traits,
#' otherwise: list with entries for each person: matrix of SEs: rows = blocks, columns = traits
#' @export
#'
info2se <- function(infos, summed=TRUE, var.out=FALSE, s.prior=NULL) {

  #SEs from infos (within a person)
  .info2se <- function(infos, summed=TRUE, var.out) {
    if (isTRUE(summed)) infos <- matrix(colSums(infos, dims=1), dim(infos)[2], dim(infos)[3])
    if (length(dim(infos))==3) {
      #several blocks
      if(isTRUE(var.out)) {
        t(apply(infos, 1, function(info) diag(MASS::ginv(info))))
      } else {
        t(apply(infos, 1, function(info) sqrt(diag(MASS::ginv(info)))))
      }
    #summed / only one block
    } else if(is.null(s.prior)) {
      if(isTRUE(var.out)) {
        diag(MASS::ginv(infos))
      } else {
        sqrt(diag(MASS::ginv(infos)))
      }
    } else {
      if(isTRUE(var.out)) {
        diag(MASS::ginv(infos + s.prior))
      } else {
        sqrt(diag(MASS::ginv(infos + s.prior)))
      }
    }
  }

  if (is.list(infos)) {
    ses <- lapply(infos, .info2se, summed=summed, var.out=var.out)
  } else {
    ses <- .info2se(infos, summed=summed, var.out=var.out)
  }

  if(isTRUE(summed)) {
    return(do.call(rbind, ses))
  } else {
    return(ses)
  }
}
