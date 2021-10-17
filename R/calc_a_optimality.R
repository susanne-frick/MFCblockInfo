#' calculate A-optimality
#'
#' @param infos block information as returned from calc.block.infos
#'
#' @return numeric vector of A-optimality for each person
#' @export
#'
#' @examples
calc.a.optimality <- function(infos) {
  rowSums(info2se(infos, var.out=T), na.rm=T)
}
