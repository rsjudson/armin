#--------------------------------------------------------------------------------------
#
# print out the receptor tree
#
#--------------------------------------------------------------------------------------
receptorTree <- function(supermatrix=SUPERMATRIX_3) {
  printCurrentFunction()
  
  n0 <- dim(supermatrix)[1]
  cat("======================================\n")
  cat(">>> Ntotal: ",n0,"\n")
  
  temp1 <- supermatrix[!is.element(supermatrix[,"maximum.receptor"],"None"),]
  n1 <- dim(temp1)[1]
  cat("======================================\n")
  cat(">>> N1 (actives):     ",n1,"\n")
  cat(">>> N0-N1 (inactive): ",(n0-n1),"\n")
  
  temp2 <- temp1[temp1[,"specificity_T"]>0.5,]
  n2 <- dim(temp2)[1]
  cat("======================================\n")
  cat(">>> N2 (high T):   ",n2,"\n")
  cat(">>> N1-N2 (low T): ",(n1-n2),"\n")
  
  temp3 <- temp2[temp2[,"specificity_Z"]>0.5,]
  n3 <- dim(temp3)[1]
  cat("======================================\n")
  cat(">>> N3 (high Z):   ",n3,"\n")
  cat(">>> N2-N3 (low Z): ",(n2-n3),"\n")
  
  cat("======================================\n")
  rec.list <- sort(unique(temp3[,"maximum.receptor"]))
  nrec <- length(rec.list)
  rcount <- vector(length=nrec,mode="integer")
  for(i in 1:nrec) {
    receptor <- rec.list[i]
    temp4 <- temp3[is.element(temp3[,"maximum.receptor"],receptor),]
    rcount[i] <- dim(temp4)[1]
  }
  ix <- sort(rcount,index.return=T,decreasing=T)$ix
  for(i in 1:nrec) {
    index <- ix[i]
    if(rcount[index]>=5) cat(rec.list[index]," \t ",rcount[index],"\n")
  }
}