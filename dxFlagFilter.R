#--------------------------------------------------------------------------------------
#
# explore the flag filtering
#
#--------------------------------------------------------------------------------------
dxFlagFilter <- function() {
  printCurrentFunction()
  
  col.name <- "maximum.receptor"
  file <- "../output/flag_filter_ztfiltered.txt"
  s <- paste("Assay\tN\tmean.in\tmean.out\tp.value\n",sep="")
  cat(s,file=file,append=F)
  cat(s)
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>0.5,]
  supermatrix <- temp[temp[,"specificity.T"]>0.5,]
  
  for(i in 1:NASSAY) {
    flag.col <- paste(ASSAY.LIST[i],"_flags",sep="")
    aname <- paste("A",i,sep="")
    mask <- supermatrix[,col.name]
    vals.in <- supermatrix[is.element(mask,aname),flag.col]
    vals.out <- supermatrix[is.element(mask,"Agonist"),flag.col]
    mean.in <- mean(vals.in)
    mean.out <- mean(vals.out)
    
    n.in <- length(vals.in)		
    if(n.in>=5) {
      result <- ks.test(vals.in,vals.out,alternative="less",exact=F)
      p.val <- result$p.value
      s <- paste(aname,"\t",n.in,"\t",format(mean.in,digits=2),"\t",format(mean.out,digits=2),"\t",format(p.val,digits=2),"\n",sep="")
      cat(s,file=file,append=T)
      cat(s)
    }
  }
}		
