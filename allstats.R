#--------------------------------------------------------------------------------------
#
# calcuate the number of chemicals in activity classes
#
#===============================================================================
# Pathway-specific - start
#===============================================================================
#--------------------------------------------------------------------------------------
allstats <- function() {
  printCurrentFunction()
  flush.console()
  
  nrange <- 3
  rangemin <- c(0.05,0.01, -1)
  rangemax <- c(  100, 0.05, 0.01)
  results <- as.data.frame(matrix(nrow=(nrange+1),ncol=1+NRECEPTOR))
  rnames <- c("R1 (Agonist)","R2 (Antagonist)")
  for(i in 3:NRECEPTOR) rnames <- c(rnames,paste("R",i,sep=""))
  rnames <- c(rnames,"R1 (Agonist) zcut","R2 (Antagonist) zcut")
  names(results) <- c("AUC Range","R1 (Agonist)","R2 (Antagonist)",rnames[3:(NRECEPTOR)])
  for(i in 1:nrange) {
    results[i,1] <- paste(rangemin[i]," to ",rangemax[i])
    for(j in 1:NRECEPTOR) {
      receptor <- paste("R",j,sep="")
      arindex <- arIndex(receptor)
      nickname <- arindex$nickname
      col <- paste("AUC.",nickname,sep="")
      temp <- SUPERMATRIX[,col]
      temp <- temp[temp>=rangemin[i]]
      temp <- temp[temp<rangemax[i]]
      results[i,j+1] <- length(temp)
    }
  }
  results[4,1]<- "Specific"
  
  i <- nrange+1
  #names(results) <- c("Range","Agonist","Antagonist","R3","R4","R5","R6","R7","A1","A2","A3","A4","A5","A7","A8","A9","A10","A11","A12","A13")
  for(j in 1:NRECEPTOR) {
    aname <- names(results)[j+1]
    temp <- SUPERMATRIX[is.element(SUPERMATRIX[,"maximum.receptor"],aname),]
    results[4,(j+1)] <- dim(temp)[1]
  }
  print(results)
  outfile <- paste("../output/allchem_ranges_",ALPHA,"_",RS0,".txt",sep="")
  write.table(results,file=outfile, row.names=F, append=FALSE, quote=F, sep = "\t")
}
