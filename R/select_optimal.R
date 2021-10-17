#' Select Optimal Block Combination
#'
#' @param info.sum summary of block information as returned from calc.info.block
#' @param traits.grid grid of traits for which block information was calculated
#' @param K integer number of blocks in the pool
#' @param K.final integer, final number of blocks
#'
#' @return list: solved should equal 0 if a solution was found, ind.opt contains the indices of selected blocks
#' @export
#'
#' @examples
select.optimal <- function(info.sum, traits.grid, info.start=NULL, K, K.final, weights.grid=NULL, constraint.list=NULL) {

  if(is.null(info.start)) info.start <- rep(0, nrow(traits.grid))
  #number of decision variables (blocks + criterion)
  M <- K + 1

  if(is.null(weights.grid)) {
    #across blocks (for weights on grid points)
    info.sum.pool <- rowSums(info.sum) + info.start
    #relative levels
    info.sum.pool <- info.sum.pool/info.sum.pool[which(rowSums(traits.grid==0)==ncol(traits.grid))]
  } else {
    info.sum.pool <- weights.grid
  }

  ####----------- create model --------------####
  lpmodel <- lpSolveAPI::make.lp(0, M)
  #specifications
  lpSolveAPI::lp.control(lpmodel, sense="max")

  #variable types
  lpSolveAPI::set.type(lpmodel, columns=1:K, type="binary")
  lpSolveAPI::set.type(lpmodel, columns=M, type="real")
  lpSolveAPI::set.bounds(lpmodel, lower=rep(0,K), upper=rep(1,K), columns=1:K)
  lpSolveAPI::set.bounds(lpmodel, lower=0, columns=M)

  #constraint on test length
  #sum of all decision variables = K.final
  lpSolveAPI::add.constraint(lpmodel, rep(1,K), "=", K.final, indices=1:K)

  #additional constraints
  if(isFALSE(is.null(constraint.list))) {
    for(l in 1:length(constraint.list$right)) lpSolveAPI::add.constraint(lpmodel, constraint.list$left[,l], constraint.list$operator[l],
                                                                        constraint.list$right[l], indices=1:K)
  }

  #relative accuracy at the grid points
  #information summed across selected blocks = M*relative information for this trait in item pool
  for (g in 1:nrow(info.sum)) {
    lpSolveAPI::add.constraint(lpmodel, c(info.sum[g,],-info.sum.pool[g]), ">=", info.start[g], indices=c(1:K,M))
  }

  #objective function
  lpSolveAPI::set.objfn(lpmodel, 1, indices = M)
  #write model to check syntax
  # lpSolveAPI::write.lp(lpmodel, "simulation/model.lp", type="lp")

  #solve equations
  reslp <- solve(lpmodel)
  if(reslp==0) { #should return 0 (optimal solution found)
    #get decision variables
    sel <- lpSolveAPI::get.variables(lpmodel)
    ind.opt <- which(sel==1)
  } else {
    ind.opt <- rep(NA, K.final)
  }
  return(list(solved=reslp, ind.opt=ind.opt))
}
