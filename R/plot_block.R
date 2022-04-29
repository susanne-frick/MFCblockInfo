#' plot block information
#'
#' @param which.blocks vector, indices for blocks to plot
#' @param info information for all blocks as returned from calc.info.block()
#' @param K integer, number of blocks
#' @param loads matrix of item loadings
#' @param trait.names names for traits, if NULL defaults to Trait 1, Trait 2 etc.
#' @param ... further parameters passed to plot3D::persp3D
#'
#' @return plot with reduction in SEs, traits in rows, fixed trait in columns
#'
#' @examples
plot.block <- function(which.blocks, info, K, loads, trait.names=NULL,
                       col="skyblue1", bty="u", col.grid="grey45",
                       r=1, phi=30, shade=.2, ltheta=60, ticktype="detailed",
                       zlab="\n\nR2", border=NA,
                       nticks=6, scale=TRUE,
                       ...){
  if(is.null(trait.names)) trait.names <- paste("Trait", 1:ncol(loads))

  nb <- nrow(loads)/K
  traits.which.blocks <- create.traits.blocks(loads, which.blocks=1:K, nb=nb)

  for(b in which.blocks) {
    #extract which traits are measured by this block
    ind.info <- which(apply(info$pairs, 1, setequal, traits.which.blocks[b,]))

    #set par: rows=traits, columns=fix.levels, 0 margins
    withr::local_par(mfrow=c(length(ind.info),
                             ifelse(nb>2, length(info$variedlevels$fix.levels), 1)),
                     mai=c(0.5,0,0,0), oma=c(0,0,2,0))

    #info about which trait
    for(tr in traits.which.blocks[b,]) {
      main.name <- paste("Info about", trait.names[tr],"for Block",b)
      #levels of fix.dim
      for(ind in ind.info) {
        #calculate r2
        r2 <- calc.info.block.r2(info.all=info$info[[ind]], wo.blocks = b)

        ####----- plot ----------####
        trait.names.varied <- trait.names[info$pairs[ind,1:3]]

        z.lim <- range(r2[,tr])
        if(nb > 2) {
          for(th3 in info$variedlevels$fix.levels){
            z.mat <- matrix(r2[info$gridnb[,3]==th3, tr],
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
           title(sub=paste(trait.names.varied[3], "=", th3), line=1, adj=.55)
          }
        } else {
          z.mat <- matrix(r2[,tr],
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
        }
      }
      mtext(main.name, side=3, outer=T, cex=1)
    }
  }
}
