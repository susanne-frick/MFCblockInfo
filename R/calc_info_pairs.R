calc.info.pairs <- function(traits, int, loads, uni, design.mat) {
  
  npairs <- K*choose(nb,2)
  
  comp.mat.a <- design.mat %*% load.mat
  comp.uni <- design.mat^2 %*% uni
  comp.load <- apply(comp.mat.a, 1, function(rw) tcrossprod(rw, rw), simplify=FALSE)
  
  all.infos <- vector("list", nrow(traits))
  for (n in 1:length(all.infos)) all.infos[[n]] <- array(dim=c(K,ncol(loads),ncol(loads)))
  
  for (n in 1:nrow(traits)) {
    inner.tirt <- (comp.mat.a%*%traits[n,] + as.vector(gamma.true))/sqrt(as.vector(design.mat^2 %*% uni))
    inner.info <- dnorm(inner.tirt)^2/(pnorm(inner.tirt)*(1-pnorm(inner.tirt)))
    
    info.pairs <- lapply(1:npairs, function(ind, cuni, cload, iinfo) 1/comp.uni[ind] * comp.load[[ind]] * inner.info[ind],
                         cuni=comp.uni, cload=comp.load, iinfo=inner.info)
    all.infos[[n]] <- array(t(matrix(unlist(info.pairs), ncol(load.mat)^2, npairs)), dim=c(npairs, ncol(load.mat), ncol(load.mat)))
  }
  
  return(all.infos)
}