#' generate item parameters
#'
#' @param design.load design matrix of questionnaire, rows=items, columns=traits, 1=positively keyed, -1=negatively keyed
#' @param K scalar, number of blocks
#' @param nb scalar, block size
#' @param load.range vector, giving the range for factor loadings
#' @param int.range vector, giving the range for item intercepts
#'
#' @export
#' @return
sim.items <- function(design.load, K, nb, load.range, int.range) {
  design.mat <- create.design.mat(K=K, nb=nb)
  
  stopifnot(
    "All traits must be measured by at least one item." =
    !any(apply(design.load, 2, function (x) all(x == 0)))
  )
  
  if(nb > 1) {
    repeat{
      loads <- runif(K*nb, min=min(load.range), max=max(load.range))
      load.mat <- loads * design.load
      if(qr(design.mat %*% load.mat)$rank==ncol(design.load)) {
        break
      }
    }
  } else  {
    loads <- runif(K*nb, min=min(load.range), max=max(load.range))
    load.mat <- loads * design.load
  }
  #draw other item parameters, combine to table items
  #table: mean, loading, uniqueness
  items <- data.frame(matrix(nrow=K*nb, ncol=3))
  colnames(items) <- c("u.mean","loads","uni")
  items$u.mean <- runif(K*nb, min=min(int.range), max=max(int.range))
  items$loads <- loads

  #uniqueness = 1-load^2
  items$uni <- 1-items$loads^2

  return(items)
}
