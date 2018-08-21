#--------------------------------------------------------------------------------------
#
# build histograms of all tops
#
#--------------------------------------------------------------------------------------
topHist <- function(to.file=F) {
  printCurrentFunction()
  
  if(to.file) {
    fname <- "../plots/top_hist.pdf"
     pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,4,1))
  for(assay in ASSAY.LIST) {
    x <- MAT.t[,assay]
    hist(x,main=assay)
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}