#' internal: create matrix of sequential comparisons
#'
#' @param nb block size
#'
#' @return matrix of sequential comparisons (last row is unimportant)
#'
#' @examples create.tr(3)
#' @source SafirsMFC
create.tr <- function(nb) {
  Tr=matrix(0,nb,nb) ## initialisiere T_n_b (0 wird in alle Zellen eingetragen)
  diag(Tr)=1 ## setze Werte in der Diagonalen gleich 1
  Tr[nb,]=1 ## setze Werte in der letzten Zeile gleich 1
  Tr[(1:(nb-1))*(nb+1)]=-1 ## setze Werte direkt oberhalb der Diagonalen auf -1
  return(Tr)
}
