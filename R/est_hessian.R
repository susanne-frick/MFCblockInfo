#' Estimate traits based on genuine likelihood
#'
#' @param FUN function to compute response probability
#' @param responses matrix of block responses, rows = persons, columns = blocks. Responses should be given as indices for rank orders, corresponding to the columns in perms.
#' @param int vector of pair intercepts (i.e., intercepts for binary outcomes of pairwise comparisons)
#' @param loads matrix of item loadings, rows = items, columns = traits
#' @param uni matrix of item uniquenesses (diagonal)
#' @param perms matrix of permutations (i.e., rank orders). Can be obtained from calling permute()
#' @param SE logical. Obtain standard errors from generalized inverse of the negative hessian at the log-likelihood? defaults to TRUE.
#' @param lh.fun function to calculate likelihood across blocks. Defaults to lh.
#' @param starts matrix of starting values for the latent traits, rows = persons, columns = traits. If NULL, all starting values are zero. Defaults to NULL.
#' @param box numeric vector of length 1. Box constraints for the latent traits are set as $\pm$ box for all traits. Defaults to 3.
#' @param ... additional arguments passed to FUN.
#'
#' @return list with 5 entries: traits = matrix of point estimates for the latent traits, row = persons, columns = traits. 
#' ses = matrix of standard errors for the trait estimates, if SE = FALSE, all entries are NA. 
#' errors, warns, messages = vectors of any errors, warnings and messages that occured during estimation, in the order of their occurence, 
#'
#' @examples
est.hessian <- function(FUN, traits, responses, int, loads, uni, perms, lh.fun=lh, ...) {

  nb <- nrow(perms)
  K <- nrow(loads)/nb
  bi <- create.block.ind(K,nb)
  bi_int <- create.blocks.int(K,nb)
  perms_int <- create.perms.int(nb, perms=perms)
  perms_order <- apply(perms, 2, order)
  Tr <- create.tr(nb)

  ses <- matrix(NA, nrow(responses), ncol(loads))

  all.infos.obs <- vector("list", nrow(traits))
  for (n in 1:length(all.infos.obs))  all.infos.obs[[n]] <- matrix(NA,ncol(loads),ncol(loads))
  
  errors <- NULL
  warns <- NULL
  messages <- NULL

  #loop over persons
  for(j in 1:nrow(responses)) {
    tryCatch({
    result <- optimHess(par=traits[j,],fn=lh.fun, lhb.fun=FUN,
                 control = list(fnscale=-1),
                 responsesj=responses[j,], loads=loads, int=int, uni=uni, bi=bi, bi_int=bi_int,
                 perms_int=perms_int, Tr=Tr, perms=perms_order, ...)
    all.infos.obs[[j]] <- -result
    ses[j,] <- sqrt(diag(MASS::ginv(-result)))
    }, error=function(e){
      errors <<- c(errors, conditionMessage(e))
      ses[j,] <- NA
    }, warning=function(w)
      warns=c(warns, conditionMessage(w)))
   }
  return(list("traits"=traits, "ses"=ses, "info" = all.infos.obs, "errors"=errors, "warns"=warns))
}
