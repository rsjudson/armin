#--------------------------------------------------------------------------------------
#
# Calculate the least squares solution for one chemical
#
#===============================================================================
# Pathway-specific - update pathway-specific line color, etc. behavior
#===============================================================================
#--------------------------------------------------------------------------------------
LSoneVa <- function(code="C80057",cname="BPA",strength="Strong",dir="../input/CRref",mode="AT",do.debug=F) {
  #TmatVa()
  file <- paste(dir,"/CRMAT_",code,".txt",sep="")
  if(mode=="AT") {
    file <- paste(dir,"/CRMAT_",code,"_",mode,".txt",sep="")
  }
  adata <- read.table(file,header=T,sep="\t")
  adata <- adata[,AMASK==1]
  
  casrn <- CHEMS[code,"casrn"]
  nconc <- length(CONCLIST)
  nassay <- sum(AMASK)
  rmat <- matrix(nrow=nconc,ncol=nassay)
  rmat[] <- 0
  adata <- as.data.frame(t(adata))
  adata <- cbind(adata,TMAT)
  conc.names <- c()
  for(i in 1:NCONC) conc.names <- c(conc.names,paste("C",i,sep=""))
  t.names <- c()
  for(i in 1:NRECEPTOR) t.names <- c(t.names,paste("T",i,sep=""))
  t.names <- t.names[RMASK==1]
  names(adata) <- c(conc.names,t.names)
  nuse <- NCONC
  
  qc.string <- paste("Chemical QC [")
  qc.temp <- as.data.frame(TOX21.QC[is.element(TOX21.QC[,"code"],code),])
  if(dim(qc.temp)[1]>0) {
    for(q in 1:dim(qc.temp)[1]) qc.string <- paste(qc.string,qc.temp[q,"grade"])
  }
  qc.string <- paste(qc.string,"]")
  qc.string=""
  # plot the interpolated assay values
  temp <- adata[1,1:NCONC]
  a1 <- as.numeric(temp[1:nuse])
  xlist <- CONCLIST
  cols <- c("black","black","black","green","green","hotpink","cyan","cyan","cyan","orange","orange","violet","violet")
  anames <- c()
  for(i in 1:NASSAY) anames <- c(anames,paste("A",i,sep=""))
  ##CHANGE?
  ltys <- c(1,2,3, 1,1,2,2,3,3, 1,2, 1,2, 1,2,1, 1,2, 1,2, 1,2)
  cols <- cols[AMASK==1]
  ltys <- ltys[AMASK==1]
  anames <- anames[AMASK==1]
  plot(a1~xlist,xlim=c(0.01,200),log="x",ylim=c(0,1),cex.lab=1.2,cex.axis=1.2,xlab="Concentration (uM)",ylab="Efficacy - fraction(PC)",lwd=3,col=cols[1],lty=ltys[1],main=paste(casrn,":",cname,"\n",qc.string),type="n")
  text(x=0.01,y=0.9,strength,cex=1,pos=4)
  counter <- 1
  
  cmed <- CYTOTOX[code,"cytotox_median_log"]
  cmad <- CYTOTOX[code,"global_mad"]
  clow <- cmed - 3*cmad
  cytotox.median <- 10**(cmed)
  cytotox.lower <- 10**(clow)
  
  if(cytotox.median<1000) {
    rect(cytotox.lower,0,1000,1,col="gray80",border="black")
    lines(c(cytotox.median,cytotox.median),c(0,1),lwd=3,col="red")
  }
  
  cautions.by.code <- CAUTION.FLAGS[is.element(CAUTION.FLAGS[,"code"],code),]
  
  for(i in 1:nassay) {
    assay <- ASSAY.LIST[i]
    max.conc <- MAT.maxconc[code,assay]
    temp <- adata[i,1:NCONC]
    ai <- as.numeric(temp[1:nuse])
    cautions.by.assay <- as.data.frame(cautions.by.code[is.element(cautions.by.code[,"assay"],assay),])
    activity <- MAT.hitcall[code,assay]
    if(!is.na(max.conc)) {
      iuse0 <- interpolateConc(xlist,ai,max.conc)$iuse
      ai <- as.numeric(temp[1:iuse0])
      xi <- xlist[1:iuse0]
    }
    else {
      xi <- xlist
    }
    if(activity==1 && dim(cautions.by.assay)[1]>0) {
      clist <- cautions.by.assay[,"flag"]
      doit <- F
      pch <- 4
      #            if(is.element("Replicate mismatch",clist)) {doit <- T; pch <- 4}
      if(is.element("Multiple points above baseline, inactive",clist)) {doit <- T; pch <- 2}
      if(is.element("Borderline inactive",clist)) {doit <- T; pch <- 2}
      if(is.element("Borderline active",clist)) {doit <- T; pch <- 4}
      if(is.element("Only one conc above baseline, active",clist)) {doit <- T; pch <- 4}
      if(is.element("Only highest conc above baseline, active",clist)) {doit <- T; pch <- 4}
      if(is.element("Gain AC50 < lowest conc & loss AC50 < mean conc",clist)) {doit <- T; pch <- 4}
      #            if(is.element("Noisy data",clist)) {doit <- T; pch <- 4}
      if(is.element("Includes potential flare region points",clist)) {doit <- T; pch <- 8}
      if(is.element("Includes potential chemical plate interlace points",clist)) {doit <- T; pch <- 8}
      
      lines(ai~xi,lwd=2,col=cols[i],lty=ltys[i])
      if(doit) {
        print(cautions.by.assay[,c("name","assay","flag")])
        flush.console()
        points(ai~xi,col="black",pch=pch,cex=1)
      }
    }
    else {
      lines(ai~xi,lwd=2,col=cols[i],lty=ltys[i])
    }
  }
  #browser()
  #------------------------------------------
  # now plot the receptor curves
  #------------------------------------------
  if(do.debug) browser()
  nreceptor <- sum(RMASK)
  resmat <- as.data.frame(matrix(nrow=NCONC,ncol=nreceptor))
  allrnames <- c()
  for(i in 1:NRECEPTOR) allrnames <- c(allrnames,paste("R",i,sep=""))
  names(resmat) <- allrnames[RMASK==1]
  resmat[] <- 0
  start <- vector(mode="numeric",length=nreceptor)
  lwr <- vector(mode="numeric",length=nreceptor)
  upr <- vector(mode="numeric",length=nreceptor)
  start[] <- 0
  lwr[] <- 0
  upr[] <- 1
  for(i in 1:nuse) {
    concname <- paste("C",i,sep="")
    A <- adata[,c(concname,t.names)]
    if(i>1) start <- res$par
    if(do.debug) res <- optim(par=start,f=AFRVa,A=A,method="L-BFGS-B",lower=lwr,upper=upr,control=list(maxit=1))
    else res <- optim(par=start,f=AFRVa,A=A,method="L-BFGS-B",lower=lwr,upper=upr,control=list(maxit=2000))        #cat(concname,"  eps/alpha: ",format(res$value/ALPHA,digits=2),"\n")
    for(j in 1:nreceptor) resmat[i,j] <- res$par[j]*SCALE
    if(res$convergence!=0) cat(i,"Convergence: ",res$convergence," Calls: ",res$counts," residual: ",res$value," : ",res$message,"\n")
    #browser()
  }
  # plot the predicted R values
  ltys <- c( 1,     1,     1,       1,       1,        2,     1,        1,        2)
  cols <- c("blue","red", "black", "green", "hotpink","cyan","orange", "violet", "violet")    
  rnames <- c()
  for(i in 1:NRECEPTOR) rnames <- c(rnames,paste("R",i,sep=""))
  cols <- cols[RMASK==1]
  rnames <- rnames[RMASK==1]
  aucr1 <- receptorScore(resmat[,1],filter=F,use.slope=F)
  aucr2 <- receptorScore(resmat[,2],filter=F,use.slope=F)
  
  ret <- pseudoAC(code,"AC50",aucr1,aucr2,F)$median.value
  ac50 <- 10**(ret)
  lines(c(ac50,ac50),c(0,1.5),col="green",lwd=3)
  
  subtitle <- paste("Agonist:",format(aucr1,digits=2)," Antagonist:",format(aucr2,digits=2))
  plot(resmat[1:nuse,1]~xlist,xlim=c(0.01,200),log="x",ylim=c(0,1.5),cex.lab=1.2,cex.axis=1.2,xlab="Concentration (uM)",ylab="Receptor Score",type="n",lwd=3,col="blue",main=paste(casrn,":",cname,"\n",subtitle))
  #cytotox.median <- CYTOTOX[code,"cytotox_median_um"]
  #cytotox.lower <- CYTOTOX[code,"cytotox_lower_bound_um"]
  if(cytotox.median<1000) {
    rect(cytotox.lower,0,1000,1.5,col="gray80",border="black")
    lines(c(cytotox.median,cytotox.median),c(0,1.5),lwd=3,col="red")
  }
  lines(c(ac50,ac50),c(0,1.5),col="green",lwd=3)
  rlist <- c(1,2,3,4,5,6,7,18,19)
  for(i in 1:NRECEPTOR0) lines(resmat[1:nuse,rlist[i]]~xlist,lwd=3,col=cols[i],lty=ltys[i])
  text(x=0.01,y=1.4,strength,cex=1,pos=4)
  
  ytop <- 0.95
  dy <- 0.1
  
  return(resmat)
}
