calc.info.pairs <- function(traits, int, loads, uni, design.mat, responses = NULL) {
  npairs <- nrow(design.mat)
  
  if(is.null(responses)) {
    responses <- matrix(1, nrow(traits), npairs)
    expected.only <- TRUE
  } else {
    expected.only <- FALSE
  }
  
  comp.mat.a <- design.mat %*% loads
  comp.uni <- design.mat^2 %*% uni
  comp.load <- lapply(1:nrow(comp.mat.a), function(rind, ca) tcrossprod(ca[rind,], ca[rind,]), ca=comp.mat.a)
  
  all.infos.exp <- all.infos.obs <- vector("list", nrow(traits))
  for (n in 1:length(all.infos.exp)) all.infos.exp[[n]] <- all.infos.obs[[n]] <- array(dim=c(npairs, ncol(loads), ncol(loads)))
  
  for (n in 1:nrow(traits)) {
    inner.tirt <- (comp.mat.a%*%traits[n,] + as.vector(int))/sqrt(as.vector(comp.uni))
    
    #expected information
    inner.info.exp <- dnorm(inner.tirt)^2/(pnorm(inner.tirt)*(1-pnorm(inner.tirt)))
    inner.info.exp[is.infinite(inner.info.exp)] <- 1e-100 #number close to zero
    #info becomes infinite when dividing by zero because prob is close to 1
    info.pairs.exp <- lapply(1:npairs, function(ind, cuni, cload, iinfo) 1/cuni[ind] * cload[[ind]] * iinfo[ind],
                         cuni=comp.uni, cload=comp.load, iinfo=inner.info.exp)
    all.infos.exp[[n]] <- array(t(matrix(unlist(info.pairs.exp), ncol(loads)^2, npairs)), dim=c(npairs, ncol(loads), ncol(loads)))
    
    #observed information
    inner.info.obs <- dnorm(inner.tirt)^2/(pnorm(inner.tirt)^responses[n,] * (1-pnorm(inner.tirt))^(1-responses[n,]) )
    inner.info.obs[is.infinite(inner.info.obs)] <- 1e-100
    info.pairs.obs <- lapply(1:npairs, function(ind, cuni, cload, iinfo) 1/cuni[ind] * cload[[ind]] * iinfo[ind],
                         cuni=comp.uni, cload=comp.load, iinfo=inner.info.obs)
    all.infos.obs[[n]] <- array(t(matrix(unlist(info.pairs.obs), ncol(loads)^2, npairs)), dim=c(npairs, ncol(loads), ncol(loads)))
  }
  
  if(isTRUE(expected.only)) {
    return(all.infos.exp)
  } else {
    return(list("expected" = all.infos.exp, "observed" = all.infos.obs))
  }
}