#' Intercept indices
#'
#' create indices for intercepts
#'
#' @param K numeric, number of blocks
#' @param nb numeric, block size
#'
#' @return matrix of intercept indices for each block
#' @export
#'
#' @examples create.blocks.int(20,3)
create.blocks.int <- function(K, nb) {
  matrix(1:(K*choose(nb, 2)),choose(nb,2),K)
}
