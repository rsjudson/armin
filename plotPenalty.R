#--------------------------------------------------------------------------------------
#
# plot the penalty term
#
#--------------------------------------------------------------------------------------
plotPenalty <- function(alpha=0.1,width=1,rs0=1.5,to.file=F) {
  x <- seq(0,3,by=0.1)
  y <- x
  for(i in 1:length(x)) y[i] <- penalty(x[i],alpha,width,rs0)
  if(to.file) {
    fname <- paste("plots/penalty_",alpha,"_",width,"_",rs0,".pdf",sep="")
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  plot(y~x,xlab="R",ylab="Penalty",cex.lab=1.5, cex.axis=1.5,type="l")
  if(to.file) dev.off()
  else browser()
}
