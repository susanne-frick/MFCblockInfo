#' Likelihood Mplus
#'
#' likelihood of rank order under Thurstonian IRT model, using pair intercepts as Mplus
#'
#' @param traits vector of traits
#' @param int vector of intercepts and minus intercepts for block b
#' @param loads matrix of factor loadings for block b
#' @param uni matrix of error (co)variances for block b
#' @param y_b index for rank order in perms for block b
#' @param perms permutations (order(rank(), decreasing=F) order(, decreasing=T))
#' @param perms_int indices of intercepts (+ or -)
#' @param Tr matrix of pairwise comparisons (only consecutive ones, last row is ignored)
#' @param nb numeric, block size
#' @param algorithm algorithm to use for numerical approximation in mvtnorm::pmvnorm, defaults to TVPACK()
#' @param loga logical: return logarithmized result? defaults to TRUE
#'
#' @return likelihood
#'
#'
lhb.mplus <- function(traits,int,loads,uni,y_b,perms,perms_int,Tr,nb,algorithm=mvtnorm::TVPACK(),loga=TRUE){
  loads=as.matrix(loads) #loads should be a matrix
  udiffs <- int[perms_int[,y_b]] + (Tr %*% (loads[perms[,y_b],]%*% traits))[-nb]
  sigdiffs2 <- (tcrossprod(Tr %*% uni[perms[,y_b],perms[,y_b]], Tr))[-nb,-nb]
  out=mvtnorm::pmvnorm(lower=rep(0,nb-1), upper=rep(Inf,nb-1), mean=udiffs, sigma=sigdiffs2, algorithm=algorithm)

  ## out: likelihood for one block and one person
  if (loga==TRUE) {
    return(ifelse(out<=0, -10^100, log(out))) #small number
  } else {
    return(ifelse(out<=0, 1e-100, out)) #number close to 0
  }
}
