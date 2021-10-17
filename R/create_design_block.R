#' generate design matrix for one block
#'
#' @param nb numeric, block size
#'
#' @return matrix: rows = pairs, columns = items
#' @export
#'
#' @examples create.design.block(3)
create.design.block <- function(nb) {
  if(nb > 1) {
    pairs.b <- combn(1:nb, 2)
    d.block <- matrix(0, nb, ncol(pairs.b))
    for(p in 1:ncol(pairs.b)) d.block[pairs.b[,p],p] <- c(1,-1)
  } else {
    d.block <- matrix(1, 1, 1)
  }
  return(t(d.block))
}
