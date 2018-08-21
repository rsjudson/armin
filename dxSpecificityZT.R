#--------------------------------------------------------------------------------------
#
# explore the specificity
#
#--------------------------------------------------------------------------------------
dxSpecificityZT <- function() {
  printCurrentFunction()
  
  file <- "../output/specificity_zt.txt"
  
  temp.agonist <- SUPERMATRIX[SUPERMATRIX[,"AUC.Agonist"]>=SPECIFIC.AUC.CUTOFF.AG,]
  temp.antagonist <- SUPERMATRIX[SUPERMATRIX[,"AUC.Antagonist"]>=SPECIFIC.AUC.CUTOFF.ANT,]
  s <- "Receptor\tN\t|T|\t|T|(Comparator)\tp(T)\t|Z|\t|Z|(Comparator)\tp(Z)\t\n"
  cat(s,file=file,append=F)
  
  for(i in 3:NRECEPTOR) {
    receptor <- paste("R",i,sep="")
    arindex <- arIndex(receptor)
    nickname <- arindex$nickname
    assay.list <- arindex$assay.list
    colname <- paste("AUC.",nickname,sep="")
    
    temp.receptor <- SUPERMATRIX[SUPERMATRIX[,colname]>=SPECIFIC.AUC.CUTOFF.ANT,]
    if(nickname=="R5"|| nickname=="R6" || nickname=="A6" || nickname=="A7" || nickname=="A8" || nickname=="A9") {
      temp.receptor <- SUPERMATRIX[SUPERMATRIX[,colname]>=SPECIFIC.AUC.CUTOFF.AG,]				
    }
    
    if(dim(temp.receptor)[1]>=10) {
      nassay <- length(assay.list)
      ac50.list <- assay.list
      T.list <- assay.list
      Z.list <- ac50.list
      for(j in 1:nassay) {
        ac50.list[j] <- paste(ac50.list[j],"_logAC50",sep="")
        T.list[j] <- paste(T.list[j],"_T",sep="")
        Z.list[j] <- paste(Z.list[j],"_Zscore",sep="")
      }
      
      cat("\n==================================\n",nickname,":",dim(temp.receptor)[1],"\n==================================\n")
      
      s <- paste(nickname,"\t",dim(temp.receptor)[1],"\t",sep="")
      
      x <- temp.receptor[,T.list]
      x <- x[temp.receptor[,ac50.list]<1000000]
      y <- temp.antagonist[,T.list]
      y <- y[temp.antagonist[,ac50.list]<1000000]	
      if(nickname=="R5"|| nickname=="R6" || nickname=="A6" || nickname=="A7" || nickname=="A8" || nickname=="A9") {
        y <- temp.agonist[,T.list]
        y <- y[temp.agonist[,ac50.list]<1000000]	
      }
      if(length(x)>0 & length(y)>0){
        t.ac50 <- t.test(x=x,y=y,alternative="less")
        cat("Emax: \t",format(mean(x),digits=2),":",format(mean(y),digits=2),":",format(t.ac50$p.value,digits=2),"\n")
        s <- paste(s,format(mean(x),digits=2),"\t",format(mean(y),digits=2),"\t",format(t.ac50$p.value,digits=2),"\t",sep="")
      }
      x <- temp.receptor[,Z.list]
      x <- x[temp.receptor[,ac50.list]<1000000]
      y <- temp.antagonist[,Z.list]
      y <- y[temp.antagonist[,ac50.list]<1000000]	
      if(nickname=="R5"|| nickname=="R6" || nickname=="A6" || nickname=="A7" || nickname=="A8" || nickname=="A9") {
        y <- temp.agonist[,Z.list]
        y <- y[temp.agonist[,ac50.list]<1000000]				
      }
      if(length(x)>0 & length(y)>0){
        z.ac50 <- t.test(x=x,y=y,alternative="less")
        cat("Z: \t",format(mean(x),digits=2),":",format(mean(y),digits=2),":",format(z.ac50$p.value,digits=2),"\n")
        s <- paste(s,format(mean(x),digits=2),"\t",format(mean(y),digits=2),"\t",format(z.ac50$p.value,digits=2),"\n",sep="")
        cat(s,file=file,append=T)
      }
      #print(sort(unique(temp.receptor[,"use_category"])))
      #cat("\n\n")
      #print(sort(unique(temp.receptor[,"structure_category"])))
      # cat("\n\n")
      # print(sort(unique(temp.receptor[,"target_gene"])))
    }
  }
}
