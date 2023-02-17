#' plot testinformation
#'
#' plots testinformation for all pairwise trait combinations
#'
#' @param info information as returned from testcalc.info.plot()
#' @param loads matrix of item loadings
#' @param trait.names vector of trait names, defaults to "Trait 1","Trait 2", ...
#' @param par.mfrow mfrow argument for par(), vector with c(rows,columns) for plot
#' @param ... further parameters passed to persp3D
#'
#' @return plot with reduction in SEs, traits in rows, fixed trait in columns
#' @export
#'
plot.testinfo <- function(info, loads=NULL, trait.names=NULL, par.mfrow=c(3,4),
                          col="skyblue1", bty="u", col.grid="grey45",
                          r=1, phi=30, shade=.2, ltheta=60, ticktype="detailed",
                          xlab="\nTrait 1", ylab="\nTrait 2", zlab="\n\nSE", border=NA,
                          nticks=6, scale=TRUE,
                          ...){
  if(is.null(trait.names)) trait.names <- paste("Trait", 1:ncol(loads))

  #save original par settings
  par.orig <- list(mfrow=par()$mfrow, mai=par()$mai, oma=par()$oma)

  #set par: rows=traits, columns=fix.levels, 0 margins
  withr::local_par(mfrow=par.mfrow, mai=c(0.2,0,0,0), oma=c(0,0,2,0))
  #if par.mfrow does not fit, determine number of plots to be left empty
  #number of frames
  if (prod(par.mfrow) > nrow(info$pairs)) { #figure is larger than plots
    left <- prod(par.mfrow) %% nrow(info$pairs)
  } else if (prod(par.mfrow) < nrow(info$pairs)) { #figure is smaller than plots
    left <- prod(par.mfrow) - nrow(info$pairs) %% prod(par.mfrow)
  }

  #info about which trait
  for(tr in 1:ncol(info$pairs)) {
    main.name <- paste("SEs for", trait.names[tr])
    #z.lim <- range(do.call(c, lapply(info$info, function(i) i[,tr])))
    #levels of fix.dim
    for(ind in 1:nrow(info$pairs)) {

      ####----- plot ----------####
      trait.names.varied <- trait.names[info$pairs[ind,1:2]]

      z.lim <- range(info$info[[ind]][,tr])
      z.mat <- matrix(info$info[[ind]][, tr],
                      ncol=length(info$variedlevels$tr.levels),nrow=length(info$variedlevels$tr.levels))
      print(plot3D::persp3D(x=info$variedlevels$tr.levels,y=info$variedlevels$tr.levels,
                            z=z.mat,
                            zlim=z.lim,
                            col=col, bty=bty, col.grid=col.grid,
                            r=r, phi=phi, shade=shade, ltheta=ltheta, ticktype=ticktype,
                            xlab=trait.names.varied[1],
                            ylab=trait.names.varied[2],
                            zlab=zlab, border=border,
                            nticks=nticks, scale=TRUE, ...))
      if(ind==1) {
        mtext(main.name, side=3, outer=T, cex=1)
      } else if ((prod(par.mfrow) < nrow(info$pairs)) & (ind %in% prod(par.mfrow) + 1)) {
        mtext(main.name, side=3, outer=T, cex=1)
      }
    }
    if (prod(par.mfrow) != nrow(info$pairs)) for(l in 1:left) plot.new()
  }
}
