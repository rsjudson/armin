#--------------------------------------------------------------------------------------
#
# explore the cross-talk in tox21
#
#--------------------------------------------------------------------------------------
dxTox21Crosstalk <- function() {
  printCurrentFunction()
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>0.5,]
  supermatrix <- temp[temp[,"specificity.T"]>0.5,]
  
  col.name <- "maximum.receptor"
  file <- "../output/tox21_crosstalk.txt"
  s <- paste("AssaySet\tReceptor\tTox1.ref.mean\tTox21.receptor.mean\tp.Tox21\n",sep="")
  cat(s,file=file,append=F)
  
  assay.set <- "BLA.Agonist"
  assays <- names(TOXCAST.TESTED)
  tox21.assays <- assays[grep("TOX21",assays)]
  tox21.assays <- tox21.assays[grep("Agonist",tox21.assays)]
  tox21.assays <- tox21.assays[grep("BLA",tox21.assays)]
  tox21.assays <- tox21.assays[!is.element(tox21.assays,"TOX21_AR_BLA_Agonist_ratio")]
  z.tox21 <- TOXCAST.ZMAT[CODE.LIST,tox21.assays]
  z.tox21[is.na(z.tox21)] <- 0
  z.tox21[z.tox21<3] <- 0
  z.tox21[z.tox21>0] <- 1
  hits.tox21 <- rowSums(z.tox21)
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.tox21)
  names(temp)[dim(temp)[2]] <- "TOX21.BLA.Agonist.nonAR.hits.HIZ"
  SUPERMATRIX <<- temp
  
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  tox21.ref <- as.numeric(hits.tox21[ref.codes])
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    tox21.receptor <- as.numeric(hits.tox21[receptor.codes])
    if(length(tox21.receptor)>=5) {
      result <- ks.test(tox21.receptor,tox21.ref,alternative="less",exact=F)
      s <- paste(assay.set,"\t",receptor,"\t",format(mean(tox21.ref),digits=2),"\t",format(mean(tox21.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
  
  assay.set <- "BLA.Antagonist"
  assays <- names(TOXCAST.TESTED)
  tox21.assays <- assays[grep("TOX21",assays)]
  tox21.assays <- tox21.assays[grep("Antagonist",tox21.assays)]
  tox21.assays <- tox21.assays[grep("BLA",tox21.assays)]
  tox21.assays <- tox21.assays[!is.element(tox21.assays,"TOX21_AR_BLA_Antagonist_ratio")]
  z.tox21 <- TOXCAST.ZMAT[CODE.LIST,tox21.assays]
  z.tox21[is.na(z.tox21)] <- 0
  z.tox21[z.tox21<3] <- 0
  z.tox21[z.tox21>0] <- 1
  hits.tox21 <- rowSums(z.tox21)
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.tox21)
  names(temp)[dim(temp)[2]] <- "TOX21.BLA.Antagonist.nonAR.hits.HIZ"
  SUPERMATRIX <<- temp
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  tox21.ref <- as.numeric(hits.tox21[ref.codes])
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    tox21.receptor <- as.numeric(hits.tox21[receptor.codes])
    if(length(tox21.receptor)>=5) {
      result <- ks.test(tox21.receptor,tox21.ref,alternative="less",exact=F)
      s <- paste(assay.set,"\t",receptor,"\t",format(mean(tox21.ref),digits=2),"\t",format(mean(tox21.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
  
  assay.set <- "LUC.Agonist"
  assays <- names(TOXCAST.TESTED)
  tox21.assays <- assays[grep("TOX21",assays)]
  tox21.assays <- tox21.assays[grep("Agonist",tox21.assays)]
  tox21.assays <- tox21.assays[grep("LUC",tox21.assays)]
  tox21.assays <- tox21.assays[!is.element(tox21.assays,"TOX21_AR_LUC_MDAKB2_Agonist")]
  z.tox21 <- TOXCAST.ZMAT[CODE.LIST,tox21.assays]
  z.tox21[is.na(z.tox21)] <- 0
  z.tox21[z.tox21<3] <- 0
  z.tox21[z.tox21>0] <- 1
  hits.tox21 <- rowSums(z.tox21)
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.tox21)
  names(temp)[dim(temp)[2]] <- "TOX21.LUC.Agonist.nonAR.hits.HIZ"
  SUPERMATRIX <<- temp
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  tox21.ref <- as.numeric(hits.tox21[ref.codes])
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    tox21.receptor <- as.numeric(hits.tox21[receptor.codes])
    if(length(tox21.receptor)>=5) {
      result <- ks.test(tox21.receptor,tox21.ref,alternative="less",exact=F)
      s <- paste(assay.set,"\t",receptor,"\t",format(mean(tox21.ref),digits=2),"\t",format(mean(tox21.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
  
  assay.set <- "LUC.Antagonist"
  assays <- names(TOXCAST.TESTED)
  tox21.assays <- assays[grep("TOX21",assays)]
  tox21.assays <- tox21.assays[grep("Antagonist",tox21.assays)]
  tox21.assays <- tox21.assays[grep("LUC",tox21.assays)]
  tox21.assays <- tox21.assays[!is.element(tox21.assays,"TOX21_AR_LUC_MDAKB2_Antagonist2")]
  z.tox21 <- TOXCAST.ZMAT[CODE.LIST,tox21.assays]
  z.tox21[is.na(z.tox21)] <- 0
  z.tox21[z.tox21<3] <- 0
  z.tox21[z.tox21>0] <- 1
  hits.tox21 <- rowSums(z.tox21)
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.tox21)
  names(temp)[dim(temp)[2]] <- "TOX21.LUC.Antagonist.nonAR.hits.HIZ"
  SUPERMATRIX <<- temp
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  tox21.ref <- as.numeric(hits.tox21[ref.codes])
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    tox21.receptor <- as.numeric(hits.tox21[receptor.codes])
    if(length(tox21.receptor)>=5) {
      result <- ks.test(tox21.receptor,tox21.ref,alternative="less",exact=F)
      s <- paste(assay.set,"\t",receptor,"\t",format(mean(tox21.ref),digits=2),"\t",format(mean(tox21.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
  outfile <- "../output/superMatrix_ATG_NVS_Tox21.csv"
  write.csv(SUPERMATRIX,file=outfile, row.names=F)
}	