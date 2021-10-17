#' Uniqueness matrix
#'
#' @param mplus.pars Mplus parameters, read in with MplusAutomation
#' @param item.identifier regular expression that is part of all item names
#'
#' @return matrix of estimated uniquenesses (diagonal)
#' @export
#'
uni2mat <- function(mplus.pars, item.identifier) {
  diag(abs(mplus.pars[grepl(".WITH",mplus.pars$paramHeader) & grepl(item.identifier,mplus.pars$param),]$est))
}
