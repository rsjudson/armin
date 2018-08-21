#--------------------------------------------------------------------------------------
#
# plot the range of AUC in the agonist channel (R1) as a function of expected potency
#
#--------------------------------------------------------------------------------------
refchemDist <- function(to.file=F,zcut=F) {
  file <- paste("../output/refchem_AUC_",ALPHA,"_",RS0,".txt",sep="")
  if(zcut) file <- paste("../output/refchem_AUC_",ALPHA,"_",RS0,"_zcut.txt",sep="")
  rdata <- read.table(file,sep="\t",header=T,stringsAsFactors=F,comment.char="",quote="\"")
  bnames <- c("Inactive","Very Weak","Weak","Moderate","Strong")
  nbins <- length(bnames)
  y <- rdata[,"R1"]
  groups <- y
  groups[] <- 0
  refout <- REFCHEMS
  refout <- cbind(refout,refout[,1])
  names(refout)[dim(refout)[2]] <- "Agonist.AUC"
  refout <- cbind(refout,refout[,1])
  names(refout)[dim(refout)[2]] <- "Antagonist.AUC"
  refout[,"Agonist.AUC"] <- 0
  refout[,"Antagonist.AUC"] <- 0
  
  for(i in 1:dim(rdata)[1]) {
    chem <- rdata[i,1]
    expect <- as.character(REFCHEMS[is.element(REFCHEMS[,1],chem),"agonist_potency"])
    if(!is.na(expect)) {
      if(expect=="Inactive") groups[i] <- 1
      if(expect=="Very Weak") groups[i] <- 2
      if(expect=="Weak") groups[i] <- 3
      if(expect=="Moderate") groups[i] <- 4
      if(expect=="Strong") groups[i] <- 5
    }
    R1 <- rdata[i,"R1"]
    R2 <- rdata[i,"R2"]
    refout[which.max(refout[,1]==chem),"Agonist.AUC"] <- R1
    refout[which.max(refout[,1]==chem),"Antagonist.AUC"] <- R2
  }
  outfile <- paste("../output/refchem_scored_",ALPHA,"_",RS0,".txt",sep="")
  if(zcut) outfile <- paste("../output/refchem_scored_",ALPHA,"_",RS0,"_zcut.txt",sep="")
  write.table(refout,file=outfile, row.names=F, append=FALSE, quote=F, sep = "\t")
  
  y <- y[groups>0]
  groups <- groups[groups>0]
  
  if(to.file) {
    fname <- paste("../plots/refchem_activity_ranges_",ALPHA,"_",RS0,".pdf",sep="")
    if(zcut) fname <- paste("../plots/refchem_activity_ranges_",ALPHA,"_",RS0,"_zcut.pdf",sep="")
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,0.1))
  main <- "Agonist Score (R1) vs. Reference Activity Class"
  if(zcut) main <- "Agonist Score (R1) vs. Reference Activity Class (Z>3)"
  browser()
  boxplot(y~groups,ylab="Agonist Score",xlab="Activity Class",names=bnames,cex.axis=1.2,cex.lab=1.2,main=main)
  if(to.file) dev.off()
  else browser()
}
