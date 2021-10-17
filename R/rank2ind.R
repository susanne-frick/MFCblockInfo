rank2ind <- function(y_b) {
  #transform y_b to index
  if (length(y_b)>1) y_b <- which(apply(perms, 2, function(perm) all(perm==y_b)))
  return(y_b)
}
