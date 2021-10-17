#' Create grid of traits
#'
#' @param tr.levels vector of trait levels
#' @param ntraits integer, number of traits
#'
#' @return dataframe with all combinations of trait levels for ntraits
#' @export
#'
#' @examples create.grid(seq(-1,1,.5), 3)
create.grid <- function(tr.levels, ntraits) {
  tr.list <- vector("list", ntraits)
  for(tr in 1:ntraits) tr.list[[tr]] <- tr.levels
  grid.nb <- expand.grid(tr.list)
  grid.nb
}
