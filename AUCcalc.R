#--------------------------------------------------------------------------------------
#
# calculate the area under the curve
#
#--------------------------------------------------------------------------------------
AUCcalc <- function(mode="ref",zcut=F) {
  file <- paste0("../output/",mode,"chem_resmat_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/",mode,"chem_resmat_",ALPHA,"_",RS0,"_zcut.xlsx")
  rdata <- read.xlsx(file)
  cmethod <- HEATMAP.CMETHOD
  nlevel <- 25
  print(dim(rdata))
  name.list <- c("code","casrn","name")
  for(i in 1:sum(RMASK)) name.list <- c(name.list,paste0("R",i))
  row <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) <- name.list
  omat <- NULL
  nmax <- dim(rdata)[1]
  for(i in 1:nmax) {
    code <- rdata[i,1]
    cname <- CHEMS[code,"name"]
    casrn <- CHEMS[code,"casrn"]
    counter <- 3
    row[] <- NA
    row[1,"code"] <- code
    row[1,"casrn"] <- casrn
    row[1,"name"] <- cname
    for(j in 1:sum(RMASK)) {
      istart <- (j-1)*NCONC + 4
      iend <- j*NCONC + 3
      x <- rdata[i,istart:iend]
      use.slope <- T
      if(j==1 || j==2) use.slope <- F
      auc <- receptorScore(x,method=2,do.print=F,use.slope=use.slope)
      row[1,3+j] <- auc
      #if(code=="C52806538") {
      #  browser()
      #  cat(j,auc,names(row)[j+3],use.slope,"\n")
      #}
    }
    omat <- rbind(omat,row)
    
  }
  file <- paste0("../output/",mode,"chem_AUC_",ALPHA,"_",RS0,".xlsx")
  if(zcut) file <- paste0("../output/",mode,"chem_AUC_",ALPHA,"_",RS0,"_zcut.xlsx")
  write.xlsx(omat,file)
}
