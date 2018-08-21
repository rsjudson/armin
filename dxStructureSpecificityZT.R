#--------------------------------------------------------------------------------------
#
# explore the specificity based on structure
#
#--------------------------------------------------------------------------------------
dxStructureSpecificityZT <- function(super=F,cutoff=0.5) {
  printCurrentFunction()
  
  col.name <- "maximum.receptor"
  scol.name <- "structure_category"
  file <- "../output/structure_category_specificity_zt.txt"
  if(super) {
    scol.name <- "structure_super_category"
    file <- "../output/structure_super_category_specificity_zt.txt"
  }
  
  temp <- SUPERMATRIX[SUPERMATRIX[,"specificity.Z"]>cutoff,]
  supermatrix <- temp[temp[,"specificity.T"]>cutoff,]
  
  txt <- TxT(1,2,3,4)
  s <- paste("Structure_category\tReceptor\t",txt$title,"\n",sep="")
  cat(s,file=file,append=F)
  str.list <- sort(unique(supermatrix[,scol.name]))
  rec.list <- sort(unique(supermatrix[,col.name]))
  nstr <- length(str.list)
  nrec <- length(rec.list)
  for(i in 2:nstr) {
    str.class <- str.list[i]
    for(j in 1:nrec) {
      receptor <- rec.list[j]
      str.mask <- supermatrix[,scol.name]
      str.mask[] <- 0
      str.mask[is.element(supermatrix[,scol.name],str.class)] <- 1
      rec.mask <- supermatrix[,col.name]
      rec.mask[] <- 0
      rec.mask[is.element(supermatrix[,col.name],receptor)] <- 1
      rec.mask[!is.element(supermatrix[,col.name],receptor)] <- -1
      str.mask <- as.numeric(str.mask)
      rec.mask <- as.numeric(rec.mask)
      
      str.mask <- str.mask[rec.mask!=0]
      rec.mask <- rec.mask[rec.mask!=0]
      rec.mask[rec.mask== -1] <- 0
      
      a <- sum(str.mask*rec.mask)
      b <- sum(str.mask*(1-rec.mask))
      c <- sum((1-str.mask)*rec.mask)
      d <- sum((1-str.mask)*(1-rec.mask))
      txt <- TxT(a,b,c,d)
      if(a>=4) {
        s <- paste(str.class,"\t",receptor,"\t",txt$sval,"\n")
        cat(s,file=file,append=T)
        cat(s)
      }
    }
  }
}