#--------------------------------------------------------------------------------------
#
# Calculate the least squares solution for all chemicals
#
#--------------------------------------------------------------------------------------
LSall <- function(to.file=F,mode="AT",zcut=F) {
  dir <- "../input/CRall"
  if(zcut) dir <- "../input/CRall.zcut"
  smask <- paste(AMASK,collapse="")
  TmatVa()
  par(mfrow=c(3,2),mar=c(4,4,4,1))
  if(to.file) {
    fname <- paste("../plots/allchem_perf_",ALPHA,"_",RS0,".pdf",sep="")
    if(zcut) fname <- paste("../plots/allchem_perf_",ALPHA,"_",RS0,"_zcut.pdf",sep="")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    par(mfrow=c(3,2),mar=c(4,4,4,1))
  }
  nreceptor <- sum(RMASK)
  name.list <- c("code","casrn","name")
  for(j in 1:nreceptor) for(k in 1:NCONC) name.list <- c(name.list,paste0("R",j,"_C",k))
  row <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) <- name.list
  omat <- NULL
  
  nmax <- dim(CHEMS)[1]
  for(i in 1:nmax) {
    iuse <- i
    code <- CHEMS[iuse,"code"]
    cname <- CHEMS[code,"name"]
    casrn <- CHEMS[code,"casrn"]
    cat(code,"[",cname,"] ",i,":",nmax,"\n")
    flush.console()
    resmat <- LSoneVa(code,cname,"",dir,mode)
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
  
  file <- paste0("../output/allchem_resmat_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/allchem_resmat_",ALPHA,"_",RS0,"_zcut.xlsx")
  write.xlsx(omat,file)
  
}
