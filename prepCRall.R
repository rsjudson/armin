#--------------------------------------------------------------------------------------
#
# prep the concentration-response matrix for all chemicals
#
#--------------------------------------------------------------------------------------
prepCRall <- function(mode="AT",zcut=F) {
  printCurrentFunction()
  dir <- "../input/CRall"
  if(zcut) dir <- "../input/CRall.zcut"
  code.list <- CHEMS[,"code"]
  for(code in code.list) {
    if(!is.na(code)) prepCR(code,dir,mode,zcut)
  }
}
