#--------------------------------------------------------------------------------------
#
# calculate the specificity.score
# NK: add specificity flag based on additional Tox21 data (MDAKB2 antagonist, high R-1881)
#
#--------------------------------------------------------------------------------------
specificityScores <- function(code,smatrix) {
  specr <- smatrix[code,"maximum.receptor"]
  if(specr=="None") return(list(spec.t=0,spec.zscore=0,spec.summary=0))
  logac50 <- smatrix[code,"pseudo.logAC50.median"]
  zmed <- CYTOTOX[code,"cytotox_median_log"]
  zmad <- CYTOTOX[code,"global_mad"]
  
  Z <- -(logac50-zmed) /zmad
  #browser()
  t.list <- NULL
  z.list <- NULL
  for(i in 1:NASSAY) {
    tname <- paste(ASSAY.LIST[i],"_T",sep="")
    tvalue <- as.numeric(smatrix[code,tname])
    if(tvalue>0) t.list <- c(t.list,tvalue)
    zname <- paste(ASSAY.LIST[i],"_Zscore",sep="")
    zvalue <- as.numeric(smatrix[code,zname])
    if(tvalue>0) z.list <- c(z.list,zvalue)
  }
  t.value <- median(t.list)
  z.value <- median(z.list)

  t.cut <- 50
  z.cut <- 3
  t.w <- 10.0
  z.w <- 0.2

  t.val <- sigmoid(t.value,t.cut,t.w)
  z.val <- sigmoid(z.value,z.cut,z.w)
  spec.score <- t.val * z.val
  
  return(list(spec.t=t.val,spec.zscore=z.val,spec.summary=spec.score))
}