#' permute
#'
#' generate all possible permutations (without replacement)
#'
#' @param v vector to permute
#'
#' @return matrix, permuted elements in columns
#' @export
#'
#' @examples permute(1:3)
#' @source https://www.r-bloggers.com/learning-r-permutations-and-combinations-with-base-r/
#'
permute <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- cbind(X, rbind(v[i], permute(v[-i])))
    X
  }
}
