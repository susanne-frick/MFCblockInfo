#' compute testinformation for plots
#'
#' computes testinformation for all pairwise trait combinations, other traits are fixed to fix.level.others
#'
#' @param tr.levels vector, trait levels that vary
#' @param fix.level.others numeric, level for traits fixed for this plot, defaults to 0
#' @param K numeric, number of blocks
#' @param loads matrix of item loadings
#' @param which.blocks vector, subset of blocks to calcuate testinfo for, defaults to 1:K
#' @param ... FUN to calculate information and further parameters passed to FUN
#'
#' @return information for relevant trait combinations for all blocks, and for which trait pairs (combinations) information is
#'
calc.testinfo.plot <- function(tr.levels, fix.level.others=0, K, loads, which.blocks, ...){

  nb <- nrow(loads)/K
  if(is.null(which.blocks)) which.blocks <- 1:K
  traits.which.blocks <- create.traits.blocks(loads, which.blocks=which.blocks, nb=nb)
  traits.test <- unique(c(traits.which.blocks))
  #order by index
  traits.test <- traits.test[order(traits.test)]
  pairs <- t(combn(traits.test, 2))

  if(length(traits.test)<4) {
    pairs.others <- cbind(pairs, apply(pairs, 1, function(p, trt) trt[! trt %in% p],trt=traits.test))
  } else {
    pairs.others <- cbind(pairs, t(apply(pairs, 1, function(p, trt) trt[! trt %in% p],trt=traits.test)))
  }

  #create grid of traits if traits are not given
  tr.list <- vector("list", 2)
  for(tr in 1:2) tr.list[[tr]] <- tr.levels
  grid.nb <- expand.grid(tr.list)

  info <- vector("list", nrow(pairs.others))
  for(po in 1:nrow(pairs.others)) {

    traits.b <- matrix(NA, nrow(grid.nb), ncol(loads))
    traits.b[,pairs.others[po,1:2]] <- as.matrix(grid.nb)[,1:2]
    if(length(traits.test) > 2) traits.b[,pairs.others[po,-c(1:2)]] <- fix.level.others

    info[[po]] <- info2se(calc.info.block(loads=loads, traits=traits.b, which.blocks=1:K, K=K, nb=nb, ...), summed=T)
  }
  return(list("info"=info, "pairs"=pairs.others, "gridnb"=grid.nb,
              "variedlevels"=list("tr.levels"=tr.levels, "fix.level.others"=fix.level.others)))
}
