#' Testinfo 1D
#'
#' calculate testinfo for one trait varied, levels for other traits are randomly sampled
#' from multivariate normal distribution with covariance sigma, averaged across sample of size n
#'
#' @param tr.levels vector of trait levels to vary on x-axis
#' @param sigma matrix of trait correlations
#' @param n numeric, random samples to draw for other traits
#' @param seed numeric, optional, random seed
#' @param which.blocks vector, blocks for which to calculate testinfo
#' @param ... FUN for block info and other arguments passed to FUN
#'
#' @return list of 1) matrix of standard errors ("ses") and 2) vector of varied levels ("tr.levels")
#' @export
#'
calc.testinfo.1d <- function(tr.levels, sigma, n=20, seed=NULL, which.blocks, ...) {

  if (is.null(seed)==FALSE) set.seed(seed)

  #sample traits
  traits.fix <- lapply(1:ncol(sigma), function(fd, trl, s, n) {
                       sample.traits(tr.levels=trl, fix.dim=fd, sigma=s, n=n)
  }, trl=tr.levels, s=sigma, n=n)

  #calculate info
  info.fix <- vector("list", length=ncol(sigma))
  for(tr in 1:ncol(sigma)) info.fix[[tr]] <- calc.info.block(traits=traits.fix[[tr]], ...)

  #testinfo
  ses.fix <- lapply(info.fix, info2se, summed=T)
  #average across levels of fix.dim
  traits.ses.fix <- matrix(NA, length(tr.levels), ncol(sigma))
  for (tr in 1:ncol(sigma)) traits.ses.fix[,tr] <- tapply(ses.fix[[tr]][,tr], traits.fix[[tr]][,tr], mean)

  return(list("ses"=traits.ses.fix, "variedlevels"=tr.levels))
}
