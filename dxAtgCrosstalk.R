#--------------------------------------------------------------------------------------
#
# explore the cross-talk in attagene
#
#--------------------------------------------------------------------------------------
dxAtgCrosstalk <- function() {
  printCurrentFunction()
  
  col.name <- "maximum.receptor"
  file <- "../output/atg_crosstalk.txt"
  assays <- ASSAY.LIST
  atg.assays <- assays[grep("ATG_",assays)]
  trans.assays <- atg.assays[grep("TRANS_",atg.assays)]
  trans.assays <- trans.assays[!is.element(trans.assays,"ATG_AR_TRANS_up")]
  
  z.trans <- MAT.z[CODE.LIST,trans.assays]
  z.trans[is.na(z.trans)] <- 0
  z.trans[z.trans<3] <- 0
  z.trans[z.trans>0] <- 1
  hits.trans <- rowSums(z.trans)
  
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.trans)
  names(temp)[dim(temp)[2]] <- "ATG.nonAR.hits.HIZ.TRANS"
  SUPERMATRIX <<- temp
  
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  trans.ref <- as.numeric(hits.trans[ref.codes])
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>0.5,]
  supermatrix <- temp[temp[,"specificity.T"]>0.5,]
  
  s <- paste("Receptor\ttrans.ref.mean\ttrans.receptor.mean\tp.trans\n",sep="")
  cat(s,file=file,append=F)
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    trans.receptor <- as.numeric(hits.trans[receptor.codes])
    result <- ks.test(trans.receptor,trans.ref,alternative="less",exact=F)
    s <- paste(s,receptor,"\t",format(mean(trans.ref),digits=2),"\t",format(mean(trans.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
    cat(s,file=file,append=T)
    cat(s)
    
  }
  outfile <- "../output/superMatrix_ATG.csv"
  write.csv(SUPERMATRIX,file=outfile, row.names=F)
}	
