#' select blocks based on greedy algorithm
#'
#' @param infos block information as returned from calc.block.info()
#' @param FUN function to summarize block information
#' @param traits.grid matrix, traits info was computed for
#' @param K integer, number of blocks in pool
#' @param K.start integer, number of blocks in initial questionnaires. It is assumed that these are the first K.start blocks.
#' @param K.final integer, number of blocks in final questionnaire
#'
#' @return numeric vector of indices of selected blocks
#' @export
#'
#' @examples
select.greedy <- function(infos, FUN, traits.grid, K, K.start, K.final, maximize=F, weights.grid=NULL, constraint.list=NULL) {
  min.FUN <- ifelse(maximize, which.max, which.min)
  if(is.null(weights.grid)) {
    #summarize information for weights on traits/grid points
    infos.all <- FUN(infos)
    if(maximize) {
      infos.all <- infos.all/infos.all[[which(rowSums(traits.grid==0)==ncol(traits.grid))]]
    } else {
      infos.all <- infos.all[[which(rowSums(traits.grid==0)==ncol(traits.grid))]]/infos.all
    }
  } else {
    infos.all <- weights.grid
  }
  #initialize decision vector
  sel <- ifelse(1:K <= K.start, 1, 0)
  if(is.null(constraint.list)) {
    while(sum(sel==1) < K.final) {
      ses.b <- vector("numeric", K)
      for(b in which(sel==0)) {
        info.b <- lapply(infos, function(ip) ip[c(which(sel==1),b),,])
        ses.b[b] <- mean(FUN(info.b) * infos.all) #weighted
      }
      sel[sel==0][min.FUN(ses.b[sel==0])] <- 1
    }
  } else {
    while(sum(sel==1) < K.final) {
      ses.b <- vector("numeric", K)
      con.b <- vector("numeric", K)
      for(b in which(sel==0)) {
        info.b <- lapply(infos, function(ip) ip[c(which(sel==1),b),,])
        ses.b[b] <- mean(FUN(info.b) * infos.all) * ifelse(maximize, -1, 1) #weighted
        con.b[b] <-   sum(do.call(c, lapply((colSums(constraint.list$left[sel==1,]) + constraint.list$left[b,]) - constraint.list$right, max, 0)))
      }
      #normalize critera
      ses.b <- ses.b/sum(ses.b)
      if(sum(con.b)>0) con.b <- con.b/sum(con.b)
      #combine criteria
      ses.con.b <- ses.b + con.b
      #select
      sel[sel==0][which.min(ses.con.b[sel==0])] <- 1
    }
  }
  ind.opt <- which(sel==1)
  return(ind.opt)
}
