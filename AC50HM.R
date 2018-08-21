library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#'
#' do the hierarchical clustering on the AUC matrix
#'
#--------------------------------------------------------------------------------------
AC50HM <- function(to.file=F,zcut=F) {
  printCurrentFunction()
  if(to.file) {
    fname <- "../plots/chem_AC50_heatmap.pdf"
    if(zcut) fname <- "../plots/chem_AC50_heatmap_zcut.pdf"
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  cmethod <- HEATMAP.CMETHOD
  dmat <- 6 - MAT.logac50
  rs <- rowSums(abs(dmat))
  #browser()
  if(zcut) {
    zmat <- ZSCORE
    zmat[zmat<3] <- 0
    zmat[zmat>0] <- 1
    dmat <- dmat*zmat
  }
  dmat <- dmat[rs>0,]
  nchem <- dim(dmat)[1]
  namelist <- c()
  #browser()
  for(i in 1:NASSAY) namelist <- c(namelist,paste("A",i,sep=""))
  names(dmat) <- namelist
  cat("getting ready to run heatmap\n")
  #browser()
  main <- paste("AC50 Heatmap:",nchem,"chemicals")
  if(zcut) main <- paste("AC50 Heatmap (Z>3):",nchem,"chemicals")
  result <- heatmap.2(as.matrix(dmat),
                      margins=c(5,5),
                      scale="none",
                      main=main,
                      xlab="",
                      ylab="",
                      cexCol=1.,
                      cexRow=0.1,
                      col=brewer.pal(9,"Reds"),
                      dendrogram="both",
                      Rowv=T,
                      Colv=T,
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      key.title="Key",
                      key.xlab="log(AC50)",
                      trace="none")
  #heatmap(as.matrix(dmat),margins=c(5,5),scale="none",labRow="",
  #        xlab="",ylab="",cexCol=1,cexRow=0.1,col=brewer.pal(9,"Reds"),
  #        hclustfun=function(x) hclust(d=dist(x),method=cmethod),
  #        main=main)
  if(to.file) dev.off()
  else browser()
}
