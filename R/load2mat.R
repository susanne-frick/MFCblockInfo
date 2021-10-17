#' load2mat
#'
#' Mplus loadings to matrix
#'
#' @param mplus.pars Mplus parameters, read in with Mplus Automation
#' @param design.load loading matrix, rows = items, columns = traits
#' @param design.mat MFC design matrix: rows = pairwise comparisons, columns = items
#'
#' @return loading matrix with Mplus estimates: rows = items, columns = traits
#'
#'
load2mat <- function(mplus.pars, design.load, design.mat) {
  load.tirt <- mplus.pars[grep("BY",mplus.pars$paramHeader),]$est
  #recode loadings for second items in pairwise comparison
  comp.mat.abs <- design.mat %*% abs(design.load)
  load.tirt <- load.tirt * c(comp.mat.abs[comp.mat.abs!=0])
  #remove duplicated loadings
  load.tirt <- load.tirt[seq(1,2*nrow(design.load),by=2)]

  #load.mat.tirt: design.load with estimated loadings
  load.mat.tirt <- design.load
  load.mat.tirt[load.mat.tirt!=0] <- load.tirt
  load.mat.tirt
}
