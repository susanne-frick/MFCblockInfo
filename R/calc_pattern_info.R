#' information for a single ranking pattern
#'
#' @param FUN response probability function
#' @param traits vector of traits for a single person
#' @param ... arguments passed to FUN
#'
#' @return matrix of ntraits x ntraits, information for this pattern
#' @export
#'
calc.pattern.info <- function(FUN=FUN, traits, ...) {
  prob <- FUN(traits=traits, loga=FALSE, ...)
  if(prob < 1e-10) prob <- 0
  info.obs <- -numDeriv::hessian(FUN, x=traits, ...)
  return(list("expected" = info.obs * prob, "observed" = info.obs))
}
