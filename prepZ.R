#--------------------------------------------------------------------------------------
#'
#' Prepare the z matrix
#' 
#' All concentrations are in log10(uM)
#'
#--------------------------------------------------------------------------------------
prepZ <- function() {
  printCurrentFunction()
  mat <- MAT.logac50
  cyto <- CYTOTOX
  zmat <- mat
  assay.list <- names(mat)
  for(assay in assay.list) {
    x <- mat[,assay]
    h <- MAT.hitcall[,assay]
    craw <- cyto[,"cytotox_median_log"]
    gmad <- cyto[,"global_mad"]
    z <- -(x-craw)/gmad
    z <- z*h
    zmat[,assay] <- z
    cat(assay,range(z),"\n")
    #plot(z~x,main=assay)
    #browser()
  }
  #browser()
  
  MAT.z <<- zmat
  zmat <- cbind(CHEMS,zmat)
  file <- "../input/AR_pathway_z.xlsx"
  write.xlsx(zmat,file)
  
}