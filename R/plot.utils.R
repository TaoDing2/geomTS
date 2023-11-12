
#' Put color bar legend on the right side of an existing plot.
#'
#' @param zlim charterers with initial index and end index
#' @param col colour bar
#'
#'
#' @export
#'
vertical.image.legend <- function(zlim, col){
  ## Get the current par information.  We restore these settings
  ## before we leave this function.
  starting.par.settings <- par(no.readonly=TRUE)
  ## Find information about the overall size of the figure and the
  ## margins included with that figure.
  mai <- par("mai")
  fin <- par("fin")
  ## Find total amount of figure space left over after the main plot
  ## is drawn. This really boils down to the area to the right of the
  ## plot, but under the top margin (which may hold a title) and above
  ## the bottom margin (which may hold a label for the x-axis).
  x.legend.fig <- c( 1.0-(mai[4]/fin[1]), 1.0 )
  y.legend.fig <- c( mai[1]/fin[2], 1.0-(mai[3]/fin[2]) )
  ## Now, give the portion of this area which can be used for the
  ## actual legend strip.  This means we need to leave a litle space
  ## between the strip and the main plot, and a bigger space for the
  ## labels to the right of the strip.  In the y direction, we can use
  ## all the space available.
  x.legend.plt <- c( x.legend.fig[1]+(0.08*(x.legend.fig[2]-x.legend.fig[1])),
                     x.legend.fig[2]-(0.6*(x.legend.fig[2]-x.legend.fig[1])) )
  y.legend.plt <- y.legend.fig
  ## Find cut/endpoints.  The lower zlim is the bottom endpoint; the
  ## upper zlim element is the lower endpoint.  Then, there are
  ## length(col)-1 cutpoints between these to get length(col)
  ## intervals.
  cut.pts <- seq(zlim[1], zlim[2], length=length(col)+1)
  ## Find the midpoint for each interval.  These are the values that
  ## will cause each color to be plotted in the color bar.
  z <- ( cut.pts[1:length(col)] + cut.pts[2:(length(col)+1)] ) / 2
  par( new=TRUE, pty="m", plt=c(x.legend.plt, y.legend.plt) )
  ## par( new=T, xpd=T, pty="m", plt=c(x.legend.plt, y.legend.plt) )
  image(x=1, y=z, z=matrix(z, nrow=1, ncol=length(col)),
        col=col, xlab="", ylab="", xaxt="n", yaxt="n")
  axis(4, mgp = c(3, 0.2, 0), las = 2, cex.axis=0.5, tcl=-0.1)
  box()
  ## Return par settings to what they were when we entered this
  ## function.  This is closely based on what is done in at the end
  ## of the image.plot() function in the "fields" package.
  mfg.settings <- par()$mfg
  par(starting.par.settings)
  par(mfg=mfg.settings, new=FALSE)
}


#' Automatically setup for grid of plots
#'
#' @param nplots  No. of plots. Maximum number is 64
#'
#' @export
#' @import graphics
auto.mfrow <- function(nplots) {
    if(nplots <= 3) par(mfrow=c(1, nplots))
    else if(nplots <= 4)  par(mfrow=c(2,2))
    else if(nplots <= 6)  par(mfrow=c(2,3))
    else if(nplots <= 9)  par(mfrow=c(3,3))
    else if(nplots <= 12) par(mfrow=c(3,4))
    else if(nplots <= 16) par(mfrow=c(4,4))
    else if(nplots <= 20) par(mfrow=c(4,5))
    else if(nplots <= 25) par(mfrow=c(5,5))
    else if(nplots <= 30) par(mfrow=c(5,6))
    else if(nplots <= 36) par(mfrow=c(6,6))
    else if(nplots <= 42) par(mfrow=c(6,7))
    else if(nplots <= 49) par(mfrow=c(7,7))
    else if(nplots <= 56) par(mfrow=c(7,8))
    else if(nplots <= 64) par(mfrow=c(8,8))
    else {
      stop("Too many plots")
    }
}
