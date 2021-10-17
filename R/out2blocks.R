#' out2blocks
#'
#' pairwise outcomes to block indices
#'
#' @param y data.frame/matrix of pairwise outcomes: rows = participants, columns = outcomes
#' @param blocks matrix of block indices, created from block.ind
#'
#' @return matrix of block indices: rows = participants, columns = blocks
#' @export
#'
#'
out2blocks <- function(y, blocks) apply(blocks, 1, function(b, df) apply(df[,b], 1, paste, collapse=""), df=y)
