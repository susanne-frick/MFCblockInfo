#' block.ind
#'
#' create indices for items
#'
#' @param K numeric, number of blocks
#' @param nb numeric, block size
#'
#' @return matrix, blocks in columns
#' @export
#'
#' @examples create.block.ind(20, 3)
#'
create.block.ind <- function(K, nb) t(matrix(seq_len(K*nb), nb, K))
