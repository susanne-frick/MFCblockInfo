#' nm
#'
#' rescale vector to add up to 1
#'
#' @param vec vector to rescale
#'
#' @return rescaled vector
#'
nm <- function(vec) exp(vec)/sum(exp(vec)) #vector to probabilities as in a nominal model, if vec are item parameters, they must fulfill a linear restriction to be estimable
