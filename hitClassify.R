#--------------------------------------------------------------------------------------
#'
#' load data required for the calculations
#'
#--------------------------------------------------------------------------------------
hitClassify <- function() {
  printCurrentFunction()
  temp <- MAT.z
  temp[is.na(temp)] <- -100
  temp[MAT.hitcall==0] <- -100
  temp[temp<3] <- 0
  temp[temp>3] <- 1
  HITS.Z.HI <<- rowSums(temp)   
  
  temp <- MAT.z
  temp[is.na(temp)] <- 100
  temp[MAT.hitcall==0] <- 100
  temp[temp<=0] <- 0.1
  temp[temp>=3] <- 0
  temp[temp>0] <- 1
  HITS.Z.LO <<- rowSums(temp)   
}