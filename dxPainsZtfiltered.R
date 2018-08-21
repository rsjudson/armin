#--------------------------------------------------------------------------------------
#
# explore the specificity based on z promiscuity
#
#--------------------------------------------------------------------------------------
dxPainsZtfiltered <- function(cutoff=0.5) {
  printCurrentFunction()
 
  col.name <- "maximum.receptor"
  file <- "../output/pains_ztfiltered.txt"
  
  
  temp <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),]
  zlo.ref <- as.numeric(temp[,"promiscuity.lo.Z"])
  zhi.ref <- as.numeric(temp[,"promiscuity.hi.Z"])
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>=cutoff,]
  supermatrix <- temp[temp[,"specificity.T"]>=cutoff,]
  
  txt <- TxT(1,2,3,4)
  s <- paste("Receptor\tReceptor.N\tHIZ.ref.mean\tHIZ.receptor.mean\tp.HIZ\tLO.ref.mean\tLOZ.receptor.mean\tp.LOZ\n",sep="")
  cat(s,file=file,append=F)
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    temp <- supermatrix[is.element(supermatrix[,col.name],receptor),]
    zlo.receptor <- as.numeric(temp[,"promiscuity.lo.Z"])
    zhi.receptor <- as.numeric(temp[,"promiscuity.hi.Z"])
    
    if(dim(temp)[1]>=5) {
      result <- ks.test(zlo.receptor,zlo.ref,alternative="less",exact=F)
      s <- paste(receptor,"\t",format(length(zlo.receptor),digits=2),"\t",format(mean(zlo.ref),digits=2),"\t",format(mean(zlo.receptor),digits=2),"\t",format(result$p.value,digits=2),"\t",sep="")
      result <- ks.test(zhi.receptor,zhi.ref,alternative="less",exact=F)
      s <- paste(s,format(mean(zhi.ref),digits=2),"\t",format(mean(zhi.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
}	
