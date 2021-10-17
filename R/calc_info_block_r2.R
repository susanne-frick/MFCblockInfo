#' information for one block in R2 metric
#'
#' r2 = 1 - var.all/var.wo
#'
#' @param info.all as calculated from info.block for all blocks
#' @param wo.blocks vector of indices of blocks for which to compute R2
#' @param s.prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#' @param FUN function to compute response probability (is info.all is not given)
#' @param ... other arguments passed to FUN
#'
#' @return matrix of reductions in SEs, rows = persons, columns = traits
#' @export
#'
calc.info.block.r2 <- function(info.all=NULL, wo.blocks=1, s.prior=NULL, FUN=NULL, ...) {
  #information for all blocks
  if(is.null(info.all)) info.all <- calc.info.block(FUN, ...)
  K <- dim(info.all[[1]])[1]

  #information without which.blocks
  info.wo <- lapply(info.all, function(info.p) info.p[(1:K)[! (1:K)%in% wo.blocks],,])
  #transform to variances
  var.all <- info2se(info.all, summed=T, var.out=T, s.prior=s.prior)
  var.wo <- info2se(info.wo, summed=T, var.out=T, s.prior=s.prior)
  #calculate reduction in SEs
  r2 <- 1 - var.all/var.wo
  #difference in posterior reliabilities
  #r2 <- fisherz2r(fisherz(1 - var.all) - fisherz(1 - var.wo))
  return(r2)
}
