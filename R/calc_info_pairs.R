calc.info.pairs <- function(traits, int, loads, uni, design.mat) {
  
  npairs <- nrow(design.mat)
  
  comp.mat.a <- design.mat %*% loads
  comp.uni <- design.mat^2 %*% uni
  comp.load <- lapply(1:nrow(comp.mat.a), function(rind, ca) tcrossprod(ca[rind,], ca[rind,]), ca=comp.mat.a)
  
  all.infos <- vector("list", nrow(traits))
  for (n in 1:length(all.infos)) all.infos[[n]] <- array(dim=c(npairs, ncol(loads), ncol(loads)))
  
  for (n in 1:nrow(traits)) {
    inner.tirt <- (comp.mat.a%*%traits[n,] + as.vector(int))/sqrt(as.vector(comp.uni))
    inner.info <- dnorm(inner.tirt)^2/(pnorm(inner.tirt)*(1-pnorm(inner.tirt)))
    inner.info[is.infinite(inner.info)] <- 1e-100 #number close to zero
    #info becomes infinite when deviding by zero because prob is close to 1
    info.pairs <- lapply(1:npairs, function(ind, cuni, cload, iinfo) 1/comp.uni[ind] * comp.load[[ind]] * inner.info[ind],
                         cuni=comp.uni, cload=comp.load, iinfo=inner.info)
    all.infos[[n]] <- array(t(matrix(unlist(info.pairs), ncol(loads)^2, npairs)), dim=c(npairs, ncol(loads), ncol(loads)))
  }
  
  return(all.infos)
}