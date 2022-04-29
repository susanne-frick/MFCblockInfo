#' likelihood across all questionnaire items
#'
#' across items for one person, maximum a posteriori with multivariate normal prior
#'
#' @param FUN function to compute response probability
#' @param traitsj vector of traits for one person
#' @param bi matrix of block indices as returned from create.block.ind()
#' @param bi_int matrix of intercep indices as returned from create.blocks.int()
#' @param int vector of item intercepts
#' @param loads matrix of item loadings, rows = items, columns = traits
#' @param uni matrix of item uniquenesses
#' @param perms permutations
#' @param m.prior vector of means used for multivariate normal prior
#' @param s.prior covariance matrix used for multivariate normal prior
#' @param ... other arguments passed to FUN
#'
#' @return list of length N, with entries (block,trait,trait)
#'
lh <- function(traitsj, lhb.fun, responsesj, bi, bi_int, int, loads, uni, m.prior, s.prior, est="MAP", ...) {
  #if (is.vector(traitsj)) traitsj <- t(matrix(traitsj))

  out <- 0
  #loop across blocks
  for (k in 1:nrow(bi)) {
    b <- bi[k,]
    b_int <- bi_int[,k]
    out <- out + lhb.fun(traits=traitsj, y_b=responsesj[k], int=c(int[b_int],-int[b_int]), loads=loads[b,], uni=uni[b,b],
                                   loga=TRUE, ...)
  }
  if(est=="MAP") {
    return(mvtnorm::dmvnorm(x=traitsj, mean=m.prior, sigma=s.prior, log=T) + out)
  } else {
    return(out)
  }
}
