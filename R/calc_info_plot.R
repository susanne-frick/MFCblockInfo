#' compute information for plots
#'
#' computes information needed for plot.block()
#'
#' @param tr.levels vector, trait levels that vary
#' @param fix.levels vector, trait levels that are fixed, only needed for block size > 2, therefore defaults to NULL
#' @param fix.level.others numeric, level for traits not measured by current block (if traits is not given), defaults to 0
#' @param K number of blocks
#' @param loads matrix of item loadings
#' @param which.blocks blocks to calculate information for, defaults to 1:K
#' @param ... FUN to calculate information and further parameters passed to FUN
#'
#' @return information for relevant trait combinations for all blocks, and for which trait pairs (combinations) information is
#'
#' @details traits measured by a block are varied pairwise, traits not in the pair are varied at fix.levels,
#' traits not measured by a block are fixed at fix.level.others
#'
#' @export
calc.info.plot <- function(tr.levels, fix.levels=NULL, fix.level.others=0, K, loads, which.blocks, ...){

  nb <- nrow(loads)/K
  if(is.null(which.blocks)) which.blocks <- 1:K
  #extract which blocks measure which traits
  traits.which.blocks <- create.traits.blocks(loads, which.blocks=which.blocks, nb=nb)
  traits.which.blocks.pasted <- apply(traits.which.blocks, 1, paste, collapse="-")
  #only unique block combinations
  blocks.unique <- apply(do.call(rbind, strsplit(unique(traits.which.blocks.pasted), split="-")), 2, as.numeric)
  #pairs in block combinations
  if(is.vector(blocks.unique)) {
    block.pairs <- combn(blocks.unique, 2)
  } else {
    block.pairs <- apply(blocks.unique, 1, combn, 2)
  }
  block.pairs <- matrix(block.pairs, ncol(block.pairs)*nrow(block.pairs)/2, 2, byrow=T)
  #order by trait index
  block.pairs <- apply(block.pairs, 1, function(rw) rw[order(rw)])
  block.pairs <- matrix(block.pairs, ncol(block.pairs)*nrow(block.pairs)/2, 2, byrow=T)
  #add missing traits (not measured by the block)
  block.pairs.others <- matrix(NA, nrow(block.pairs), nb)
  block.pairs.others[,1:2] <- block.pairs
  n.blocks.unique <- ifelse(is.vector(blocks.unique), 1, nrow(blocks.unique))
  if(nb > 2) {
    for(b in 1:n.blocks.unique) {
      block.pairs.b <- block.pairs[(nb*b-nb+1):(nb*b),]
      if(is.vector(blocks.unique)) {
        bub <- blocks.unique
      } else {
        bub <- blocks.unique[b,]
      }
      block.pairs.others[(nb*b-nb+1):(nb*b),3:nb] <- apply(block.pairs.b, 1, function(bpb, bu) bu[! bu %in% bpb],bu=bub)
    }
  }

  #create grid of traits
  tr.list <- vector("list", nb)
  for(tr in 1:nb) tr.list[[tr]] <- tr.levels
  if(nb > 2) tr.list[[(1:nb)[-c(1:2)]]] <- fix.levels
  grid.nb <- expand.grid(tr.list)

  info <- vector("list", nrow(block.pairs.others))
  for(bpo in 1:nrow(block.pairs.others)) {

    traits.b <- matrix(fix.level.others, nrow(grid.nb), ncol(loads))
    traits.b[,block.pairs.others[bpo,1:2]] <- as.matrix(grid.nb)[,1:2]
    traits.b[,block.pairs.others[bpo,-c(1:2)]] <- as.matrix(grid.nb)[,-c(1:2)]

    info[[bpo]] <- calc.info.block(loads=loads, traits=traits.b, which.blocks=1:K, K=K, ...)
  }
  return(list("info"=info, "pairs"=block.pairs.others, "gridnb"=grid.nb,
              "variedlevels"=list("tr.levels"=tr.levels, "fix.levels"=fix.levels,"fix.level.others"=fix.level.others)))
}
