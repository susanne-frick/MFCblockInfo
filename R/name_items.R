#' name.items
#'
#' create itemnames for a design load matrix
#'
#' @param design.load matrix: rows = items, cols = traits
#' @param prefix prefix for item names, before 01,02...10
#'
#' @return item names
#' @export
#'
#'
name.items <- function(design.load, prefix) ifelse((1:nrow(design.load))<10,paste0(prefix,"0",1:nrow(design.load)),paste0(prefix,1:nrow(design.load)))
