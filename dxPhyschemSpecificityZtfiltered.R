#--------------------------------------------------------------------------------------
#
# explore the specificity based on physchem properties
#
#--------------------------------------------------------------------------------------
dxPhyschemSpecificityZtfiltered <- function(cutoff=0.5) {
  printCurrentFunction()
  
  col.name <- "maximum.receptor"
  file <- "../output/physchem_specificity_ztfiltered.txt"
  
  filename <- "../input/ToxCast_physchem_QP_Chembl_electrophil_DFT.csv"
  physchem <- read.csv(file=filename,stringsAsFactors=F)
  
  temp <- physchem[,"species"]
  temp2 <- temp
  temp2[] <- 0
  temp2 <- as.numeric(temp2)
  temp2[] <- NA
  temp2[is.element(temp,"NEUTRAL")] <- 0
  temp2[is.element(temp,"ACID")] <- 1
  temp2[is.element(temp,"BASE")] <- 1
  physchem <- cbind(physchem,temp2)
  names(physchem)[dim(physchem)[2]] <- "Charged"
  
  
  s <- paste("Variable\tReceptor\tN.in\tN.out\tnorm.in\tnorm.out\tp.value\n")
  cat(s,file=file,append=F)
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>=cutoff,]
  supermatrix <- temp[temp[,"specificity.T"]>=cutoff,]
  #PHYSCHEM <<- physchem
  rownames(physchem) <- physchem[,"CODE"]
  #physchem <- physchem[row.names(supermatrix),]
  pnames <- names(physchem)[10:dim(physchem)[2]]
  nparam <- length(pnames)
  pclass <- pnames
  pclass[] <- "numeric"
  pclass[38] <- "character"
  pclass[39] <- "character"
  
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    rec.mask <- supermatrix[,col.name]
    rec.mask[] <- 0
    rec.mask[is.element(supermatrix[,col.name],receptor)] <- 1
    codes.in <- supermatrix[rec.mask==1,"CODE"]
    codes.out <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
    for(j in 1:nparam) {
      if(pclass[j]=="numeric" && length(codes.in)>=5) {
        param <- pnames[j]
        y.in <- physchem[codes.in,param]
        y.out <- physchem[codes.out,param]
        y.in <- y.in[!is.na(y.in)]
        y.out <- y.out[!is.na(y.out)]
        
        n.in <- length(y.in)
        n.out <- length(y.out)
        mean.in <- mean(y.in)
        mean.out <- mean(y.out)
        p.val <- 1
        if(!is.na(mean.in)) {
          if(!is.na(mean.out)) {
            if(n.in>=5 && n.out>=5) {
              res <- t.test(y.in,y.out)
              cat(param," : ",receptor,"\n")
              print(res)
              p.val <- res$p.value
            }
          }
        }
        if(n.in>=5) {
          s <- paste(param,"\t",receptor,"\t",n.in,"\t",n.out,"\t",format(mean.in,digits=2),"\t",format(mean.out,digits=2),"\t",format(p.val,digits=2),"\n",sep="")
          cat(s,file=file,append=T)
          cat(s)		
        }
      }
    }
  }
  physchem <- physchem[CODE.LIST,10:dim(physchem)[2]]
  temp <- cbind(SUPERMATRIX,physchem)
  SUPERMATRIX <<- temp
  outfile <- "../output/superMatrix_ATG_NVS_Tox21_physchem.csv"
  write.csv(SUPERMATRIX,file=outfile, row.names=F)
}		
