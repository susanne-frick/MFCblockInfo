#' MFC design matrix
#'
#' @param K numeric, number of blocks
#' @param db matrix, optional. block design as returned from create.design.block(), if NULL: computed from nb
#' @param nb numeric, optional. Block size
#'
#' @return design matrix of mfc: rows=pairwise comparisons, columns=items
#' @export
#'
#' @examples
create.design.mat <- function(K, db=NULL, nb=NULL) {
  if(is.null(db) & is.null(nb)) stop("Either db or nb should be provided.")
  if(is.null(db)) db <- create.design.block(nb)
  #number of pairs per block
  npair <- nrow(db)
  #design matrix: rows=blocks*pairs, cols=items
  design.mat <- matrix(0, K*npair, K*nb)
  #fill blocks in design.mat with block design
  for (k in 1:K) {
    design.mat[(npair*k-(npair-1)):(npair*k),(nb*k-(nb-1)):(nb*k)] <- db
  }
  return(design.mat)
}
