#' Calculate D-optimality
#'
#' @param infos block information as returned from calc.info.block()
#'
#' @return numeric vector of D-optimality for each person
#' @export
#'
#' @examples
calc.d.optimality <- function(infos) {
  do.call(c, lapply(infos, function(ip) det(colSums(ip, dims=1))))
}
