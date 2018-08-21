#--------------------------------------------------------------------------------------
#
# Calculate the least squares solution for all reference chemicals
#
#--------------------------------------------------------------------------------------
LSref <- function(to.file=F,mode="AT",zcut=F) {
  dir <- ".././input/CRref"
  if(zcut) dir <- "../input/CRref.zcut"
  smask <- paste(AMASK,collapse="")
  TmatVa()
  par(mfrow=c(3,2),mar=c(4,4,4,1))
  if(to.file) {
    fname <- paste("../plots/refchem_perf_",ALPHA,"_",RS0,".pdf",sep="")
    if(zcut) fname <- paste("../plots/refchem_perf_",ALPHA,"_",RS0,"_zcut.pdf",sep="")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    par(mfrow=c(3,2),mar=c(4,4,4,1))
  }

  nreceptor <- sum(RMASK)

  name.list <- c("code","casrn","name")
  for(j in 1:nreceptor) for(k in 1:NCONC) name.list <- c(name.list,paste0("R",j,"_C",k))
  row <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) <- name.list
  omat <- NULL
  nmax <- dim(REFCHEMS)[1]
  for(i in 1:nmax) {
    iuse <- i
    code <- REFCHEMS[iuse,"code"]
    cname <- CHEMS[code,"name"]
    casrn <- CHEMS[code,"casrn"]
    cat(code,"[",cname,"] ",i,":",nmax,"\n")
    flush.console()
    potency.string <- ""
    if(!is.na(REFCHEMS[iuse,4])) {
      potency.string <- paste("Agonist: ",REFCHEMS[iuse,4],sep="")
      if(!is.na(REFCHEMS[iuse,5])) potency.string <- paste(potency.string,"\nAntagonist: ",REFCHEMS[iuse,5],sep="")
    }
    else {
      potency.string <- paste("Antagonist: ",REFCHEMS[iuse,5],sep="")
    }
    resmat <- LSoneVa(code,cname,potency.string,dir=dir,mode)
    s <- c(code,casrn,cname)
    for(j in 1:nreceptor) {
      for(k in 1:NCONC) s <- c(s,resmat[k,j])
    }
    row[] <- s
    omat <- rbind(omat,row)
    if(!to.file) browser()
  }
  if(to.file) dev.off()
  for(i in 4:dim(omat)[2]) {
    x <- as.numeric(omat[,i])
    omat[,i] <- x
  }
  file <- paste0("../output/refchem_resmat_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/refchem_resmat_",ALPHA,"_",RS0,"_zcut.xlsx")
  write.xlsx(omat,file)
}
