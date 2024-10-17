#' simulate responses under Thurstonian IRT model
#'
#' @param traits matrix of trait levels, rows = persons, columns = traits
#' @param items matrix of item intercepts, loadings, and uniquenesses, as from sim.items()
#' @param design.load design matrix of questionnaire, rows=items, columns=traits, 1=positively keyed, -1=negatively keyed
#' @param K scalar, number of blocks
#' @param nb scalar, block size
#' @param return.index logical, return index for ranks as in permute() (TRUE) or ranks (FALSE)? defaults to TRUE. 
#' The ranks (and their indices) are interpreted in the following way: A rank of 1 means that the item in the respective position was ranked first (i.e. highest).
#'
#' @return matrix of ranks, rows = persons, columns = items
#' @export
#'
sim.responses <- function(traits, items, design.load, K, nb, return.index=TRUE) {
  #sample error from N(0,uniqueness)
  #error: matrix: rows=persons, cols=items
  error <- matrix(nrow=nrow(traits), ncol=nrow(items))
  for (i in 1:nrow(items)) {
    error[,i] <- rnorm(n=nrow(traits), mean=0, sd=sqrt(items[i,"uni"]))
  }

  load.mat <- items$loads * design.load
  ##compute utilities
  util <- t(items[,"u.mean"] + load.mat %*% t(traits)) + error
  ##utility differences
  util.diff <- util %*% t(create.design.mat(K=K, nb=nb))
  #outcomes
  y <- ifelse(util.diff > 0, 1, 0)
  # compute rank-order
  if (nb > 1) {
  y_rank <- t((apply(util, 1, function(p, bs) apply(bs, 1, function(b) order(order(p[b], decreasing=T))),bs=create.block.ind(K, nb))))
  # rank-order to index in perms
  y_rank.i <- apply(create.block.ind(K, nb), 1, function(bi, resp, pms) {
    apply(resp, 1, function(resp, b, p) which(apply(p, 2, function(perm) all(perm==resp[b]))), b=bi, p=pms)
  }, resp=y_rank, pms=permute(1:nb))
  } else {
    y_rank <- NULL
    y_rank.i <- NULL
  }
  if(isTRUE(return.index)) {
    return(y_rank.i)
  } else {
    return(list("outcomes"=y, "ranks"=y_rank, "rankindices"=y_rank.i))
  }

}
