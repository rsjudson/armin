#--------------------------------------------------------------------------------------
#
# Error function for contraint problem
#
#--------------------------------------------------------------------------------------
AFRVa <- function(x,A) {
  Ameas <- A[,1]
  nreceptor <- sum(RMASK)
  nassay <- sum(AMASK)
  F <- as.matrix(A[,2:(nreceptor+1)])
  R <- matrix(nrow=nreceptor,ncol=1)
  R[] <- x
  Apred <- F%*%R
  w <- vector(mode="numeric",length=NASSAY)
  w[] <- 1
  w <- w[AMASK==1]
  eret <- 0
  bot <- 0
  top <- 0
  mask <- Ameas
  mask[] <- 1
  mask[is.na(Ameas)] <- 0
  for(i in 1:nassay) {
    if(mask[i]==1) {
      top <- top + w[i] * (Apred[i]-Ameas[i])**2
      bot <- bot + w[i]
    }
  }
  eret <- top/bot/sum(mask)
  pen <- penalty(sum(R),ALPHA,WIDTH,RS0)
  eret <- eret+pen
  return(eret)
}
