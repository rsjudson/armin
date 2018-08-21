#--------------------------------------------------------------------------------------
#
# add the specificity score to the supermatrix
#
#--------------------------------------------------------------------------------------
addSpecificityScore <- function(smatrix=SUPERMATRIX_0) {
  printCurrentFunction()
  temp <- smatrix
  temp <- cbind(temp,temp[,dim(temp)[2]])
  names(temp)[dim(temp)[2]] <- "specificity_T"
  temp <- cbind(temp,temp[,dim(temp)[2]])
  names(temp)[dim(temp)[2]] <- "specificity_Z"

  temp[,"specificity_T"] <- 0
  temp[,"specificity_Z"] <- 0
  temp[,"specificity_score"] <- 0
  
  nchem <- dim(temp)[1]
  for(i in 1:nchem) {
    code <- temp[i,"code"]
    res <- specificityScores(code,smatrix)
    temp[i,"specificity_T"] <- res$spec.t
    temp[i,"specificity_Z"] <- res$spec.zscore
    temp[i,"specificity_score"] <- res$spec.summary
  }
  SUPERMATRIX_1 <<- temp
  file <- "../output/superMatrix_1.xlsx"
  write.xlsx(SUPERMATRIX_1,file)
  
}