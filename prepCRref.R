#--------------------------------------------------------------------------------------
#
# prep the concentration-response matrix for the reference chemicals
#
#--------------------------------------------------------------------------------------
prepCRref <- function(mode="AT",zcut=F) {
  printCurrentFunction()
  dir <- "../input/CRref"
  if(zcut) dir <- "../input/CRref.zcut"
  code.list <- REFCHEMS[,"code"]
  for(code in code.list) {
    prepCR(code,dir,mode,zcut,do.debug=F)
  }
}
