#' empirical reliability based on factor scores
#'
#' @param scores matrix of factor scores, rows = persons, columns = traits
#' @param ses matrix of SEs, same as scores
#'
#' @return vector of empirical reliabilities for each trait
#' @export
#'
calc.rel.emp <- function(scores, ses) {
  if(is.null(dim(scores))==TRUE){
    var(scores)/(var(scores) + mean(ses^2, na.rm = TRUE))
  }
  else{
  svar <- apply(scores, 2, var)
  svar/(svar + colMeans(ses^2, na.rm=T))
  }
}
