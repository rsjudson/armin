#--------------------------------------------------------------------------------------
#
# prep the concentration-response matrix for a single chemical
#
#--------------------------------------------------------------------------------------
prepCR <- function(code="C85687",dir="input/CRref.zcut",mode="AT",zcut,do.debug=F) {
  cat(code,"\n")
  flush.console()
  
  ac50 <- vector(length=NASSAY,mode="numeric")
  top <- vector(length=NASSAY,mode="numeric")
  w <- vector(length=NASSAY,mode="numeric")
  z <- vector(length=NASSAY,mode="numeric")
  emax <- vector(length=NASSAY,mode="numeric")
  maxconc <- vector(length=NASSAY,mode="numeric")
  
  for(i in 1:NASSAY) {
    assay <- ASSAY.LIST[i]
    logac50 <- MAT.logac50[code,assay]
    ac50[i] <- 10**(logac50)
    top[i] <-  MAT.t[code,assay]
    w[i] <-  MAT.w[code,assay]
    emax[i] <-  MAT.emax[code,assay]
    maxconc[i] <-  MAT.maxconc[code,assay]
    z[i] <- MAT.z[code,assay]
  }
  top[is.na(top)] <- 0
  emax[is.na(emax)] <- 0
  w[is.na(w)] <- 1
  if(do.debug) {
    print(ac50)
    print(top)
    print(w)
    print(emax)
    print(maxconc)
    print(z)
    browser()
  }
  z[z<3] <- 0
  z[z>3] <- 1
  top[top>1000] <- 0
  w[w>1000] <- 1
  top <- top/100
  top[top>1] <- 1
  cr.mat <- matrix(nrow=length(CONCLIST),ncol=NASSAY)
  cr.mat[] <- 0
  for(i in 1:length(CONCLIST)) {
    conc <- CONCLIST[i]
    for(j in 1:NASSAY) {
      ac50j <- as.numeric(ac50[j])
      tj <- as.numeric(top[j])
      wj <- as.numeric(w[j])
      emaxj <- as.numeric(emax[j])
      zj <- as.numeric(z[j])
      if(is.na(emaxj)) emaxj <- 0
      if(zj==0 && zcut) {
        ac50j <- 1000000
        tj <- 0
        if(do.debug)cat("flatten assay: ",j,"\n")
      }
      
      cr.mat[i,j] <- tj*(conc**wj/(conc**wj+ac50j**wj))
      if(ac50[j]==2000000) cr.mat[i,j] <- NA
    }
  }
  
  fname <- paste(dir,"/CRMAT_",code,"_",mode,".txt",sep="")
  write.table(cr.mat,fname,sep="\t",row.names=F)
  return(cr.mat)
}
