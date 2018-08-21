#--------------------------------------------------------------------------------------
#
# do the hierarchical clustering on the R matrix
#
#===============================================================================
# Pathway-specific - update pathway-specific assay behavior
#===============================================================================
#--------------------------------------------------------------------------------------
RHM <- function(mode="ref",to.file=F,zcut=F) {
  if(to.file) {
    fname <- paste("../plots/",mode,"chem_heatmap_",ALPHA,"_",RS0,".pdf",sep="")
    if(zcut) fname <- paste("../plots/",mode,"chem_heatmap_",ALPHA,"_",RS0,"_zcut.pdf",sep="")
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  file <- paste0("../output/",mode,"chem_resmat_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/",mode,"chem_resmat_",ALPHA,"_",RS0,"_zcut.xlsx")
  cmethod <- HEATMAP.CMETHOD
  nlevel <- 25
  print(file)
  #browser()
  rdata <- read.xlsx(file)
  
  nchem <- dim(rdata)[1]
  nconc <- 13
  #dmat <- matrix(nrow=nchem,ncol=sum(RMASK)*(nconc-1))
  counter <- 0
  
  dmat <- as.matrix(rdata[,4:dim(rdata)[2]])
  #dmat <- dmat[,1:(NRECEPTOR0*NCONC)]
  if(mode=="ref") {
    color.list <- vector(length=dim(REFCHEMS)[1],mode="character")
    color.list <- REFCHEMS[,"agonist_potency"]
    #color.list[] <- "white"
    color.list[is.element(color.list,"Strong")] <- "red"
    color.list[is.element(color.list,"Active")] <- "red"
    color.list[is.element(color.list,"Strong/Moderate")] <- "red"
    color.list[is.element(color.list,"Weak")] <- "yellow"
    color.list[is.element(color.list,"Very Weak")] <- "green"
    color.list[is.element(color.list,"Very weak")] <- "green"
    color.list[is.element(color.list,"Moderate")] <- "orange"
    color.list[is.element(color.list,"Moderate/Weak")] <- "orange"
    color.list[is.element(color.list,"Inactive")] <- "gray"
    color.list[is.element(color.list,"Negative")] <- "gray"
    color.list[is.element(color.list,"Antagonist")] <- "red"
    color.list[is.element(color.list,"Antagonist Inactive")] <- "white"
    color.list[is.element(color.list,"Cytostatic")] <- "gray"
    color.list[is.element(color.list,"SARM")] <- "yellow"
    color.list[is.element(color.list,"Agonist")] <- "green"
  }
  main <- "Normalized Response by Concentration"
  if(zcut) main <- "Normalized Response by Concentration Z>3"
  
  if(mode=="ref") {
    heatmap(as.matrix(dmat),margins=c(5,10),scale="none",labRow=rdata[,2],Colv=NA,revC=T,
            xlab="",ylab="",cexCol=0.1,cexRow=1,col=brewer.pal(9,"Reds"),
            hclustfun=function(x) hclust(d=dist(x),method=cmethod),RowSideColors=color.list,main=main)
  }
  else {
    heatmap(as.matrix(dmat),margins=c(5,10),scale="none",labRow=rdata[,2],Colv=NA,revC=T,
            xlab="",ylab="",cexCol=0.1,cexRow=1,col=brewer.pal(9,"Reds"),
            hclustfun=function(x) hclust(d=dist(x),method=cmethod),main=main)
  }
  if(to.file) dev.off()
  cl <- hclust(d=dist(dmat),method=cmethod)
  hcut <- 1
  clcut <- cutree(cl,h=hcut)
  clout <- cbind(clcut,clcut)
  clout <- cbind(clout,clcut)
  clout<- as.data.frame(clout)
  for(i in 1:length(clcut)) {
    clout[i,1] <- rdata[i,1]
    clout[i,2] <- rdata[i,2]
  }
  names(clout) <- c("code","name","Level_1")
  
  cat("Finished preping clusters for hcut: ",hcut,"\n")
  flush.console()
  
  for(hcut in 2:nlevel) {
    clcut <- cutree(cl,h=hcut)
    clout <- cbind(clout,clcut)
    names(clout)[dim(clout)[2]] <- paste("Level_",hcut,sep="")
    cat("Finished preping clusters for hcut: ",hcut,"\n")
    flush.console()
  }
  outfile <- paste("../output/",mode,"chem_clusters_",ALPHA,"_",RS0,".txt",sep="")
  if(zcut) outfile <- paste("../output/",mode,"chem_clusters_",ALPHA,"_",RS0,"_zcut.txt",sep="")
  write.table(clout,file=outfile, row.names=F, append=FALSE, quote=F, sep = "\t")
}
