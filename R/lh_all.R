#' round2
#'
#' likelihood of rank order for all items (and one person) under Thurstonian IRT model
#'
#' @param traits vector of traits
#' @param int vector of intercepts
#' @param loads matrix of factor loadings
#' @param uni matrix of error covariances
#' @param y_b rank order
#' @param perms permutations
#' @param perms_int indices of intercepts (+ or -)
#' @param lhb.fun function to compute likelihood
#'
#' @return likelihood of all
#'
lh.all <- function(th,int,loads,uni,y_b,blocks,blocks_int,perms,perms_int,lhb.fun=lhb.fun,...){
  ##### (gibt log-Likelihood eines Forced-Choice-Tests aus) ########################
  if (is.matrix(blocks)){blocks=as.data.frame(blocks)} ## R-intern: stelle sicher das blocks eine Liste ist
  if (is.matrix(blocks_int)){blocks_int=as.data.frame(blocks_int)} ## R-intern: stelle sicher das blocks eine Liste ist
  if (is.matrix(y_b)){y_b=as.data.frame(y_b)} ## R-intern: stelle sicher das blocks eine Liste ist
  out=vector("numeric",length=length(blocks)) ## initialisiere log-Likelihood des Forced-Choice Tests
  for (i in 1:length(blocks)){ ## bestimme die permsumme der log-Likelihoods über alle Itemblöcke
    sel=blocks[[i]] ## bestimme b_i
    sel_int <- blocks_int[[i]]
    out[i]=lhb.fun(traits=traits,int=int[sel_int],loads=loads[sel,],uni=uni[sel,sel],y_b=y_b[[i]],perms=perms,perms_int=perms_int,...)
  }
  return(out)
}
