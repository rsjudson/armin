#--------------------------------------------------------------------------------------
#
# Calculate the EAS quantity for each chemical
#
#--------------------------------------------------------------------------------------
EAScalc <- function(mode="FP",receptor="R2",cutoff1=0.1,cutoff2=0,eas.cutoff1=4,eas.cutoff2=6) {
  cat("EAScalc: ",receptor,"\n")
  flush.console()
  if(mode=="FP") {
    if(!exists("DMAT.FP")) {
      file <- "../input/structure_input/ToxCast_Tanimoto_matrix_REDUCED_2013_03_05.txt"
      temp <- read.table(file,sep="\t",header=T,check.names=F,quote="",comment.char="",stringsAsFactors=F)
      DMAT.FP <<- temp
      cat("DMAT.FP",dim(DMAT.FP),"\n")    	
    }
    dmat <- DMAT.FP
  }
  if(mode=="CT") dmat <- DMAT.CT
  
  if(!exists("AUC")) {
    file <- paste0("../output/allchem_AUC_",ALPHA,"_",RS0,".xlsx")
    temp <- read.xlsx(file)
    AUC <<- temp
  }
  
  mask <- as.data.frame(AUC[,receptor])
  maskA <- AUC[AUC[,receptor]>=cutoff1,"code"]
  maskB <- AUC[AUC[,receptor]<=cutoff2,"code"]
  nchem <- dim(AUC)[1]
  score <- as.data.frame(matrix(nrow=nchem,ncol=8+NRECEPTOR))
  rnames <- c()
  for(i in 1:NRECEPTOR) rnames <- c(rnames,paste("R",i,sep=""))
  names(score) <- c("code","name",rnames,"p.ttest","p.wilcox","EAS.score","EAS.class","KNN.score","KNN.class")
  for(i in 1:nchem) {
    code <- AUC[i,"code"]
    cname <- AUC[i,"name"]
    val <- AUC[i,receptor]
    if(val<1e-4) val <- 0
    score[i,"code"] <- code
    score[i,"name"] <- cname
    for(j in 1:NRECEPTOR) {
      column <- paste("R",j,sep="")
      score[i,column] <- AUC[i,column]
    }
    score[i,"p.ttest"] <- -1
    score[i,"p.wilcox"] <- -1
    score[i,"EAS.score"] <- 0
    score[i,"EAS.class"] <- ""
    score[i,"KNN.class"] <- ""
    tmaskA <- maskA[!is.element(maskA,code)]
    tmaskB <- maskB[!is.element(maskB,code)]
    tmaskA <- tmaskA[is.element(tmaskA,row.names(dmat))]
    tmaskB <- tmaskB[is.element(tmaskB,row.names(dmat))]
    
    tdA <- dmat[tmaskA,code]
    tdB <- dmat[tmaskB,code]
    if(sum(tdA)*sum(tdB)>0) {
      pttest <- t.test(tdA,tdB,alternative="greater")$p.value
      pwilco <- wilcox.test(tdA,tdB,alternative="greater")$p.value
      eas.score <- -log10(pwilco)
      score[i,"p.ttest"] <- pttest
      score[i,"p.wilcox"] <- pwilco
      score[i,"EAS.score"] <- eas.score
      if(eas.score<=eas.cutoff1 && val>=cutoff1) score[i,"EAS.class"] <- "EAS low AUC hi"
      else if(eas.score>=eas.cutoff2 && val<=cutoff2) score[i,"EAS.class"] <- "EAS hi AUC low"
      else if(eas.score>=eas.cutoff2 && val>=cutoff2 && val<cutoff1) score[i,"EAS.class"] <- "EAS hi AUC med"
    }
    score[i,"KNN.score"] <- NA
    if(is.element(code,names(dmat))) {
      dnames <- rownames(dmat)
      di <- dmat[,code]
      index <- sort(di,decreasing=T,index.return=T)$ix
      vals <- sort(di,decreasing=T,index.return=T)$x
      index <- index[1:6]
      vals <- vals[1:6]
      codesj <- dnames[index]
      vals <- vals[codesj!=code]
      codesj <- codesj[codesj!=code]
      codesj <- codesj[vals>0.7]
      if(length(codesj)>0) {
        rsj <- AUC[is.element(AUC[,"code"],codesj),receptor]
        if(length(rsj)>0) {
          knn <- max(rsj,na.rm=T)
          if(knn<1e-3) knn <- 0
          score[i,"KNN.score"] <- knn
          if(knn<=cutoff2 && val>=cutoff1) score[i,"KNN.class"] <- "KNN low AUC hi"
          else if(knn>=cutoff1 && val<=cutoff2) score[i,"KNN.class"] <- "KNN hi AUC low"
          else if(knn>=cutoff1 && val>=cutoff2 && val<cutoff1) score[i,"KNN.class"] <- "KNN hi AUC med"
        }
      }
    }
  }
  if(receptor=="R1") EAS.R1 <<- score
  if(receptor=="R2") EAS.R2 <<- score
  
  file <- paste0("../output/EAS_",mode,"_",receptor,"_",cutoff1,"_",cutoff2,"_",ALPHA,"_",RS0,".xlsx")
  write.xlsx(score,file)
}
