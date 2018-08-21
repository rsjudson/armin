#--------------------------------------------------------------------------------------
#
# receptor score or AUC
#
#--------------------------------------------------------------------------------------
receptorScore <- function(x,method=2,do.print=F,filter=F,cutoff=NA,use.slope=T) {
  nuse <- length(x)
  nuse.0 <- nuse
  if(filter==T && !is.na(cutoff)) {
    cutoff <- 10**(-cutoff) *1000000
    nuse <- which.min(CONCLIST<cutoff)
    if(do.print)cat("nuse: ",nuse,"\n")
    print("receptor.score: cutoff method no longer used")
    exit()
  }
  if(nuse<2) return(0)
  if(method==1) {
    wt <- seq(1,nuse,by=1)
    score <- sum(x[1:nuse]*wt)/sum(wt)
  }
  if(method==2) {
    if(do.print) {
      cat("==========================================\n")
      print(x)
    }
    score <- as.numeric(x[1])
    for(i in 2:nuse) {
      slope.sign <- 1
      delta <- x[i]-x[i-1]
      if(delta < -0.01) slope.sign <- -1
      if(!use.slope) slope.sign <- 1
      score <- score + slope.sign*x[i]
      if(do.print) cat("   ",i,"delta: ",format(delta,digits=3),"sign: ",slope.sign," x:",format(x[i],digits=3),"\n")
    }
    score <- score/nuse.0
    if(score<0) score <- 0
  }
  if(do.print) cat("Score: ",format(score,digits=3),"\n")
  #browser()
  return(as.numeric(score))
}
