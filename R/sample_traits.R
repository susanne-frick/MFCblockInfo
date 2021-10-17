#' sample traits from multivariate normal distribution with one trait fixed
#'
#' @param tr.levels vector: levels for fixed trait
#' @param fix.dim scalar: dimension to fix
#' @param sigma matrix: trait correlations
#' @param n scalar: number of samples for each trait level, defaults to 10
#'
#' @return traits correlating as specified in sigma, with dimension fix.dim fixed to tr.levels
#'
#' @examples sample.traits(seq(-2,2,1), fix.dim=2, sigma=diag(3), n=100)
sample.traits <- function(tr.levels, fix.dim, sigma, n=10) {
  #re-order sigma so that fix.dim is the 1st dimension

  #pairwise trait comparisons
  trait.names <- 1:ncol(sigma)
  pairs <- apply(combn(trait.names, 2), 2, paste, collapse="")
  #change traits 1 and fix.dim
  pairs.reordered <- combn(c(fix.dim,trait.names[-c(fix.dim)]), 2)
  #trait with smaller index appears on top
  pairs.reordered <- apply(pairs.reordered, 2, function(cl) cl[order(cl)])
  pairs.reordered <- apply(pairs.reordered, 2, paste, collapse="")
  cor.order <- do.call(c, lapply(pairs.reordered, function(pre, po) which(pre==po), po=pairs))

  sigma.ordered <- diag(1, ncol(sigma), ncol(sigma))
  sigma.ordered[lower.tri(sigma.ordered)] <- sigma[lower.tri(sigma)][cor.order]
  sigma.ordered[upper.tri(sigma.ordered)] <- t(sigma.ordered)[upper.tri(sigma.ordered)]

  #compute Cholesky decomposition of sigma.ordered
  R <- chol(sigma.ordered, pivot = TRUE) #pivoting allows for semi-definite matrices, but makes re-ordering necessary
  R <- R[, order(attr(R, "pivot"))]

  #randomly sample values for other traits (assuming mean=0, sd=1)
  n.tr <- length(tr.levels)
  mat <- cbind(rep(scale(tr.levels), each=n), matrix(rnorm(n.tr*n*(ncol(sigma)-1)), nrow = n.tr*n, byrow = T))

  #multiply with Cholesky to obtain specified trait correlations
  traits <- mat %*% R

  #bring fix.dim back to original position
  traits <- traits[,order(c(fix.dim,trait.names[-c(fix.dim)]))]

  #re-scale fixed trait
  traits[,fix.dim] <- traits[,fix.dim]*sd(tr.levels)

  return(traits)
}
