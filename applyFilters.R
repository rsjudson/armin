#--------------------------------------------------------------------------------------
#'
#' Apply data filters
#'
#--------------------------------------------------------------------------------------
applyFilters <- function() {
  printCurrentFunction()

  cat("Novascreen 50% x low AC50 filter\n")
  assay.list <- c("NVS_NR_hAR","NVS_NR_rAR","NVS_NR_cAR")
  for(i in 1:length(assay.list)) {
    assay <- assay.list[i]
    mask <- MAT.hitcall[,assay]
    emax <- MAT.emax[,assay]
    top <- MAT.t[,assay]
    logac50 <- MAT.logac50[,assay]
    mask[top<50] <- 0

    MAT.hitcall[,assay] <- mask
    top[mask==0] <- 0
    MAT.t[,assay] <- top
    logac50[mask==0] <- 6
    MAT.logac50[,assay] <- logac50
  }
  
  MAT.logac50 <<- MAT.logac50
  MAT.t <<- MAT.t
  MAT.hitcall <<- MAT.hitcall
}
