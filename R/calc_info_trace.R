#' Trace of information matrix
#'
#' calculates the traits (= sum of diagonal elements) of the information matrix
#'
#' @param infos block information as returned from calc.info.block()
#' @param avg logical, average result across persons? defaults to FALSE
#'
#' @return matrix of traces of information matrices, rows=persons, columns=blocks
#' @export
#'
#' @details This is the criterion of A-optimality used in computerized adaptive testing (CAT).
#' It can be used to summarize information for a block in one number, weighting all precision for all traits equally.
#' It is inversely proportional to standard errors.
#'
calc.info.trace <- function(infos, avg=F) {
  info.trace <- do.call(rbind, lapply(infos, function(ip) apply(ip, 1, function(i) sum(diag(i)))))
  if(avg==TRUE) {
    return(colSums(info.trace))
  } else {
    return(info.trace)
  }
}
