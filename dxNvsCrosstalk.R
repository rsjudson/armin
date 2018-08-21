#--------------------------------------------------------------------------------------
#
# explore the cross-talk in novascreen
#
#--------------------------------------------------------------------------------------
dxNvsCrosstalk <- function() {
  printCurrentFunction()
  
  col.name <- "maximum.receptor"
  file <- "../output/nvs_crosstalk.txt"
  assays <- names(TOXCAST.TESTED)
  nvs.assays <- assays[grep("NVS_NR",assays)]
  nvs.assays <- nvs.assays[!is.element(nvs.assays,"NVS_NR_hAR")]
  nvs.assays <- nvs.assays[!is.element(nvs.assays,"NVS_NR_cAR")]
  nvs.assays <- nvs.assays[!is.element(nvs.assays,"NVS_NR_rAR")]
  
  #	"Tox21_ERa_BLA_Agonist_ratio"	
  #	"Tox21_ERa_LUC_BG1_Agonist"	
  #	"Tox21_ERa_BLA_Antagonist_ratio"
  #	"Tox21_ERa_LUC_BG1_Antagonist"
  
  z.nvs <- TOXCAST.ZMAT[CODE.LIST,nvs.assays]
  z.nvs[is.na(z.nvs)] <- 0
  z.nvs[z.nvs<3] <- 0
  z.nvs[z.nvs>0] <- 1
  hits.nvs <- rowSums(z.nvs)
  temp <- SUPERMATRIX
  temp <- cbind(temp,hits.nvs)
  names(temp)[dim(temp)[2]] <- "NVS.nonAR.hits.HIZ"
  SUPERMATRIX <<- temp
  
  ref.codes <- SUPERMATRIX[is.element(SUPERMATRIX[,col.name],"None"),"CODE"]
  nvs.ref <- as.numeric(hits.nvs[ref.codes])
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>0.5,]
  supermatrix <- temp[temp[,"specificity.T"]>0.5,]
  
  s <- paste("Receptor\tnvs.ref.mean\tnvs.receptor.mean\tp.nvs\n",sep="")
  cat(s,file=file,append=F)
  rec.list <- sort(unique(supermatrix[,col.name]))
  nrec <- length(rec.list)
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    receptor.codes <- supermatrix[is.element(supermatrix[,col.name],receptor),"CODE"]
    nvs.receptor <- as.numeric(hits.nvs[receptor.codes])
    if(length(nvs.receptor)>=5) {
      result <- ks.test(nvs.receptor,nvs.ref,alternative="less",exact=F)
      s <- paste(receptor,"\t",format(mean(nvs.ref),digits=2),"\t",format(mean(nvs.receptor),digits=2),"\t",format(result$p.value,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
  outfile <- "../output/superMatrix_ATG_NVS.csv"
  write.csv(SUPERMATRIX,file=outfile, row.names=F)
}