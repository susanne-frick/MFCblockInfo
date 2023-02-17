#' information for one block in R2 metric
#'
#' r2 = 1 - var.all/var.wo
#'
#' @param info.all as calculated from info.block for all blocks
#' @param wo.blocks vector of indices of blocks for which to compute R2
#' @param prior matrix of prior covariances, if given, posterior information is computed, defaults to NULL
#'
#' @return matrix of reductions in SEs, rows = persons, columns = traits
#' @export
#'
calc.info.block.r2 <- function(info.all=NULL, wo.blocks=1, prior=NULL, 
                               avg = FALSE, rep = 20,
                               K.final = NULL, constraint.list = NULL, ...) {
  K <- dim(info.all[[1]])[1]

  .calc.info.block.r2 <- function(info.all, wo.blocks, prior, ...) {
    K <- dim(info.all[[1]])[1]
    #information without which.blocks
    info.wo <- lapply(info.all, function(info.p) info.p[(1:K)[! (1:K)%in% wo.blocks],,])
    #transform to variances
    var.all <- info2se(info.all, summed=T, var.out=T, prior=prior)
    var.wo <- info2se(info.wo, summed=T, var.out=T, prior=prior)
    #calculate reduction in SEs
    r2 <- 1 - var.all/var.wo
    return(r2)
  }

  r2.random.subset <- function(info.all, wo.blocks, prior, K, K.final, ...) {
    # draw rep random reference sets (fullfilling the constraints)
    random.info <- matrix(runif(K, 0, 1), 1, K)
    
    #add to constraint.list that wo.blocks should definitely be selected
    constraint.list$left <- cbind(constraint.list$left, ifelse(1:K %in% wo.blocks, 1, 0))
    constraint.list$operator <- c(constraint.list$operator, "=")
    constraint.list$right <- c(constraint.list$right, length(wo.blocks))
    
    ind.rand <- select.optimal(info.sum=random.info, traits.grid=matrix(0,1,1), K=K, K.final=K.final,
                               constraint.list=constraint.list)$ind.opt
    
    # reduce info.all to subset
    info.sub <- lapply(info.all, function(info.p, ir) info.p[ir,,], ir = ind.rand)
    # calculate r2
    r2 <- .calc.info.block.r2(info.all = info.sub, wo.blocks = which(ind.rand %in% wo.blocks), prior = prior, ...)
    return(r2)
  }
    
  if(isTRUE(avg)) {
    # calculate r2 for rep random subsets
    r2s <- replicate(rep, 
                     r2.random.subset(info.all = info.all, wo.blocks = wo.blocks, prior = prior, 
                                      K = K, K.final = K.final))
    # average
    mean.r2 <- rowMeans(r2s, dims = 2)
    # variance
    sd.r2 <- apply(r2s, c(1,2), sd)
    return(list("mean.r2" = mean.r2, "mean.sd.r2" = mean.r2 / sd.r2))
  
    } else {
    r2 <- .calc.info.block.r2(info.all = info.all, wo.blocks = wo.blocks, prior = prior, ...)
    return(r2)
  }
  
  
}
