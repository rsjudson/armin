#--------------------------------------------------------------------------------------
#
# prepare the supermatrix
#
#--------------------------------------------------------------------------------------
prepSupermatrix <- function() {
  printCurrentFunction()
  flush.console()
  file <- paste("../output/allchem_AUC_",ALPHA,"_",RS0,".xlsx",sep="")
  rdata <<- read.xlsx(file)
  result <- CHEMS
  rownames(result) <- result[,"code"]
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "EDSP.Universe"
  result[,"EDSP.Universe"] <- as.integer(is.element(CODE.LIST,EDSP.UNIVERSE[,"CODE"]))

  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "specificity_score"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "antagonist_confirmation_score"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "antagonist_confidence_score"
  
  result[,"specificity_score"] <- 0
  result[,"antagonist_confirmation_score"] <- 0
  result[,"antagonist_confidence_score"] <- 0

  result <- cbind(result,CYTOTOX[,"cytotox_median_um"])
  names(result)[dim(result)[2]] <- "cytotox_median_um"
  result <- cbind(result,CYTOTOX[,"cytotox_lower_bound_um"])
  names(result)[dim(result)[2]] <- "cytotox_lower_bound_um"
  result <- cbind(result,CYTOTOX[,"nhit"])
  names(result)[dim(result)[2]] <- "cytotox.assays.hit"
  
  for(i in 1:NRECEPTOR) {
    receptor <- paste("R",i,sep="")
    result <- cbind(result,rdata[,receptor])
    arindex <- arIndex(receptor)
    nickname <- arindex$nickname
    nickname <- paste("AUC.",nickname,sep="")
    names(result)[dim(result)[2]] <- nickname
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.logac50[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_logAC50",sep="")
  }

  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.emax[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_Emax",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.t[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_T",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.w[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_W",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.z[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_Zscore",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.maxconc[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_maxConc",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.logac10[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_logAC10",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.logacc[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_logACC",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.logacb[,assay])
    names(result)[dim(result)[2]] <- paste(assay,"_logACB",sep="")
  }
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    result <- cbind(result,MAT.logacb[,assay])
    flag.name <- paste(assay,"_flags",sep="")
    names(result)[dim(result)[2]] <- flag.name
    result[,flag.name] <- 0
    temp <- CAUTION.FLAGS[is.element(CAUTION.FLAGS[,"assay"],assay),]
    for(j in 1:NCHEM) {
      code <- CODE.LIST[j]
      temp2 <- temp[is.element(temp[,"code"],code),]
      result[code,flag.name] <- dim(temp2)[1]
    }
  }
  
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logAC50.median"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logAC50.min"
  
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logAC10.median"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logAC10.min"
  
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logACC.median"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logACC.min"
  
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logACB.median"
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "pseudo.logACB.min"
  
  result <- cbind(result,result[,dim(result)[2]])
  names(result)[dim(result)[2]] <- "maximum.receptor"
  
  ntop <- 4 + NRECEPTOR - 1
  result[,"maximum.receptor"] <- "None"
  temp <- rdata[,4:ntop]
  rs <- rowSums(temp)
  rs[rs<SPECIFIC.AUC.CUTOFF.ANT] <- 0
  rs[rs>0] <- 1
  for(i in 1:NCHEM) {
    code <- CODE.LIST[i]
    if(rs[i]==1) {		# sum of all receptor values has to be at least SPECIFIC.AUC.CUTOFF
      temp <- rdata[i,4:ntop]
      imax <- as.integer(which.max(temp))
      vimax <- temp[imax]
      if(vimax>=SPECIFIC.AUC.CUTOFF.ANT) { #  one receptor must reach this value
        
        rtemp <- paste("R",imax,sep="")
        arindex <- arIndex(rtemp)
        receptor <- arindex$nickname
        
        #temp[imax] <- 0
        #for(j in 1:length(temp)) {
        #  if(temp[j]>=SPECIFIC.AUC.CUTOFF.ANT) {
        #    rj <- paste("R",j,sep="")
        #    arindex <- arIndex(rj)
        #    nickname <- arindex$nickname
        #    if(nickname=="Agonist" || receptor=="Agonist") receptor <- "Agonist"
        #    else if(nickname=="Antagonist" || receptor=="Antagonist") receptor <- "Antagonist"
        #    #else receptor <- paste(receptor,nickname,sep=".")
        #  }
        #}
        result[i,"maximum.receptor"] <- receptor
      }
    }
  }
  # all of these are log(M) values
  for(i in 1:NCHEM) {
    code <- CODE.LIST[i]
    auc.R1 <- rdata[i,"R1"]
    auc.R2 <- rdata[i,"R2"]
    temp <- pseudoAC(code,"AC50",auc.R1,auc.R2)
    ac50 <- temp$median.value
    ac50.min <- temp$min.value
    
    temp <- pseudoAC(code,"AC10",auc.R1,auc.R2)
    ac10 <- temp$median.value
    ac10.min <- temp$min.value
    
    temp <- pseudoAC(code,"ACC",auc.R1,auc.R2)
    acc <- temp$median.value
    acc.min <- temp$min.value
    
    temp <- pseudoAC(code,"ACB",auc.R1,auc.R2)
    acb <- temp$median.value
    acb.min <- temp$min.value
    
    #cat(ac50,ac10,acc,acb,"\n")
    #if(auc.R1>0 || auc.R2>0) browser()
    if(!is.na(ac50)) {
      if(is.na(ac50.min)) ac50.min <- ac50 - 1.5
      if(ac50<1000000) {
        if(is.na(ac10)) {
          ac10 <- ac50 - 1.5
          ac10.min <- ac50.min - 1.5
        }
        if(is.na(acc)) {
          acc <- ac50 - 1.5
          acc.min <- ac50.min - 1.5
        }
        if(is.na(acb)) {
          acb <- ac50 - 1.5
          acb.min <- ac50.min - 1.5
        }
        if(ac10==6) {
          ac10 <- ac50 - 1.5
          ac10.min <- ac50.min - 1.5
        }
        if(acc==6) {
          acc <- ac50 - 1.5
          acc.min <- ac50.min - 1.5
        }
        if(acb==6) {
          acb <- ac50 - 1.5
          acb.min <- ac50.min - 1.5
        }
      }
    }
    result[i,"pseudo.logAC50.median"] <- ac50
    result[i,"pseudo.logAC10.median"] <- ac10
    result[i,"pseudo.logACC.median"] <- acc
    result[i,"pseudo.logACB.median"] <- acb
    result[i,"pseudo.logAC50.min"] <- ac50.min
    result[i,"pseudo.logAC10.min"] <- ac10.min
    result[i,"pseudo.logACC.min"] <- acc.min
    result[i,"pseudo.logACB.min"] <- acb.min
  }
  #browser()
  
  #  	rownames(EAS.R1) <- EAS.R1[,"CODE"]
  #  	rownames(EAS.R2) <- EAS.R2[,"CODE"]
  # 
  #     for(i in 1:NCHEM) {
  #         code <- CODE.LIST[i]
  #         result[i,"EAS.R1"] <- EAS.R1[code,"EAS.score"]
  #         result[i,"EAS.R1.Class"] <- EAS.R1[code,"EAS.class"]
  #         result[i,"EAS.R2"] <- EAS.R2[code,"EAS.score"]
  #         result[i,"EAS.R2.Class"] <- EAS.R2[code,"EAS.class"]
  # 	}
  
  #result <- cbind(result,NTESTED)
  #names(result)[dim(result)[2]] <- "N.assays.tested"
  result <- cbind(result,rowSums(MAT.hitcall))
  names(result)[dim(result)[2]] <- "N.assays.hit"
  result <- cbind(result,HITS.Z.HI)
  names(result)[dim(result)[2]] <- "N.assays.hit.hi.Z"
  result <- cbind(result,HITS.Z.LO)
  names(result)[dim(result)[2]] <- "N.assays.hit.lo.Z"
  
  x <- HITS.Z.HI / NASSAY
  y <- HITS.Z.LO / NASSAY
  result <- cbind(result,x)
  names(result)[dim(result)[2]] <- "promiscuity.hi.Z"
  result <- cbind(result,y)
  names(result)[dim(result)[2]] <- "promiscuity.lo.Z"
  
  SUPERMATRIX_0 <<- result
  file <- "../output/superMatrix_0.xlsx"
  write.xlsx(SUPERMATRIX_0,file)
}
