#--------------------------------------------------------------------------------------
#
# receptor score or AUC
#
#--------------------------------------------------------------------------------------
interpolateConc <- function(x,y,x0) {
  nx <- length(x)
  if(x0>=x[nx]) {
    y0 <- y[nx]
    iuse <- nx
  }
  else {
    iuse <- which.max(x0<x)
    if(iuse==0) y0 <- 0
    else if(iuse>=nx) y0 <- y[nx]
    else {
      y0 <- y[iuse]+(y[iuse+1]-y[iuse])*(log10(x0)-log10(x[iuse])) / (log10(x[iuse+1])-log10(x[iuse]))
    }
  }
  return(list(iuse=iuse,y0=y0))
}
