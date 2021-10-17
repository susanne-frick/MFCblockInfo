#' Extract traits from loading matrix
#'
#' @param loads matrix of loadings
#' @param which.blocks vector of block indices to select
#' @param nb numeric, block size
#'
#' @return matrix, rows=blocks, columns=traits
#' @export
#'
create.traits.blocks <- function(loads, which.blocks, nb) {
  K <- nrow(loads)/nb
  traits.which.blocks <- matrix(apply(loads[create.block.ind(K,nb)[which.blocks,],], 1, function(rw) which(rw!=0)), length(which.blocks),  nb)
  traits.which.blocks <- t(apply(traits.which.blocks, 1, function(rw) rw[order(rw)]))
  return(traits.which.blocks)
}
