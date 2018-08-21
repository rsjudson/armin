#--------------------------------------------------------------------------------------
#
# do the hierarchical clustering on the AUC matrix
#
#===============================================================================
# Pathway-specific - update pathway-specific assay behavior
#===============================================================================
#--------------------------------------------------------------------------------------
AUCHM <- function(mode="ref",to.file=F,zcut=F,dcut=5) {
  if(to.file) {
    fname <- paste("../plots/",mode,"chem_AUC_heatmap_",ALPHA,"_",RS0,".pdf",sep="")
    if(zcut) fname <- paste("../plots/",mode,"chem_AUC_heatmap_",ALPHA,"_",RS0,"_zcut.pdf",sep="")
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  file <- paste0("../output/",mode,"chem_AUC_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/",mode,"chem_AUC_",ALPHA,"_",RS0,"_zcut.xlsx")
  cmethod <- HEATMAP.CMETHOD
  nlevel <- 25
  
  rdata <- read.xlsx(file)
  rownames(rdata) <- rdata[,"code"]
  if(mode=="ref") rdata <- rdata[REFCHEMS[,"code"],]
  nchem <- dim(rdata)[1]
  dmat <- matrix(nrow=nchem,ncol=sum(RMASK))
  counter <- 0
  dmat <- rdata[,4:dim(rdata)[2]]
  if(mode=="ref") dmat <- dmat[,1:NRECEPTOR]
  if(mode=="all") dmat <- dmat[,1:NRECEPTOR]
  dmat.d <- dmat
  #dmat.d[dmat.d<1e-4] <- 1e-4
  dmat.d[dmat<0.05] <- 0
  dmat.d[dmat.d>0] <- 1
  rs <- rowSums(dmat.d)
  if(mode=="all") dmat <- dmat[rs>0,]
  cat("getting ready to run heatmap\n")
  flush.console()
  
  lab.col <- names(dmat)
  for(i in 1:length(lab.col)) {
    lab.col[i] <- arIndex(lab.col[i])$nickname
  }
  
  #lab.col[1] <- "Agonist"
  #lab.col[2] <- "Antagonist"
  dmat.plot <- dmat*100
  dmat.plot[dmat.plot<=dcut] <- 1
  dmat.plot <- log10(dmat.plot)
  main <- "log(AUC)"
  if(zcut) main <- "log(AUC) Z>3"
  if(mode=="ref") {
    color.list <- vector(length=dim(REFCHEMS)[1],mode="character")
    color.list[] <- "white"
    #CHANGE
    color.list[is.element(REFCHEMS[,"antagonist_potency"],c("Moderate","Moderate/Weak","Strong","Strong/Moderate","Very Weak","Weak"))] <- "red"
    color.list[is.element(REFCHEMS[,"agonist_potency"],c("Strong","Weak","Moderate"))] <- "blue"
    color.list[is.element(REFCHEMS[,"agonist_potency"],"Negative")] <- "white"

    heatmap(as.matrix(dmat.plot),margins=c(5,12),scale="none",labRow=rdata[,3],labCol=lab.col,Colv=NA,revC=T,
            xlab="",ylab="",RowSideColors=color.list,cexCol=1.2,cexRow=0.8,col=brewer.pal(9,"Reds"),main=main,
            hclustfun=function(x) hclust(d=dist(x),method=cmethod),verbose=T)
  }
  else {
    for(i in 1:NRECEPTOR) {
      receptor <- paste("R",i,sep="")
      lab.col[i] <- arIndex(receptor)$nickname
    }
    nc <- dim(dmat.plot)[1]
    main <- paste("log(AUC) chems:",nc)
    if(zcut) main <- paste("log(AUC) Z>3 chems:",nc)
    
    heatmap(as.matrix(dmat.plot),margins=c(5,5),scale="none",labRow="",labCol=lab.col,Colv=NA,revC=T,
            xlab="",ylab="",cexCol=1,cexRow=0.1,col=brewer.pal(9,"Reds"),main=main,
            hclustfun=function(x) hclust(d=dist(x),method=cmethod))
  }
  if(to.file) dev.off()
  cat("heatmap done\n")
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
  names(clout) <- c("CODE","Name","Level_1")
  
  cat("Finished prepping clusters for hcut: ",hcut,"\n")
  flush.console()
  
  for(hcut in 2:nlevel) {
    clcut <- cutree(cl,h=hcut)
    clout <- cbind(clout,clcut)
    names(clout)[dim(clout)[2]] <- paste("Level_",hcut,sep="")
    cat("Finished prepping clusters for hcut: ",hcut,"\n")
    flush.console()
  }
  outfile <- paste("../output/",mode,"chem_AUC_clusters_",ALPHA,"_",RS0,".txt",sep="")
  if(zcut) outfile <- paste("../output/",mode,"chem_AUC_clusters_",ALPHA,"_",RS0,"_zcut.txt",sep="")
  write.table(clout,file=outfile, row.names=F, append=FALSE, quote=F, sep = "\t")
  
  if(!to.file) browser()
}
