#' pair intercept indices for each permutation
#'
#' @param nb numeric, block size
#' @param perms permutations as returned from permute(). computed internally if NULL. 
#' The premutations are interpreted as ranks: A rank of 1 means that the item in the respective position was ranked first (i.e. highest).
#'
#' @return matrix: rows = intercepts, columns = permutations
#' @export
#'
#' @examples create.perms.int(4)
create.perms.int <- function(nb, perms=NULL) {
  if(nb > 1) {
  cc <- cbind(combn(1:nb,2),combn(1:nb,2)[2:1,])
  if(is.null(perms)) perms <- permute(1:nb)
  perms_int <- matrix(NA, nb-1, ncol(perms))
  for (yi in 1:ncol(perms)) {
    y_b <- order(perms[,yi], decreasing=F)
    # y_b <- perms[,yi]
    perms_int[,yi] <- do.call(c, lapply(1:(nb-1), function(ind, yb, ccc) which(apply(ccc, 2, function(perm) all(perm==yb[ind:(ind+1)]))), yb=y_b, ccc=cc))
  }
  return(perms_int)
  } else {
    return(matrix(1, 1, 1))
  }
}
