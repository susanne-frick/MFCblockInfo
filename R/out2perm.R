#' Outcomes to indices for permutations
#'
#' @param y matrix of outcomes, rows = participants, columns = pairwise outcomes
#' @param blocks matrix, created with block.ind
#' @param perms permutations
#' @param design.block block design
#'
#' @return matrix of permutation indices, rows = participants, columns = blocks
#' @export
#'
out2perm <- function(y, blocks, perms, design.block) {
  #paste outcomes within blocks
  b.y <- out2blocks(y, blocks)
  #recode patterns including NA (e.g.1NANA) as NA
  b.y[grep("NA",b.y)] <- NA

  #recodes for pattern indices
  ranks.out <- apply(perms, 2, function(perm, db) ifelse(design.block %*% perm > 0, 1, 0), db=design.block)
  ranks.out <- apply(ranks.out, 2, paste0, collapse="")

  y.ind <- apply(b.y, c(1,2), function(y, ro) ifelse(is.na(y),NA,which(y==ro)), ro=ranks.out)
  return(y.ind)
}
