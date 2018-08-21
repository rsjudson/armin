#--------------------------------------------------------------------------------------
#
# add the antagonist shift
#
#--------------------------------------------------------------------------------------
addAntagonistScore <- function(smatrix=SUPERMATRIX_2) {
  printCurrentFunction()
  temp <- smatrix
  nchem <- dim(temp)[1]
  for(i in 1:nchem) {
    code <- temp[i,"code"]
    score <- 0
    if(temp[i,"AUC.Antagonist"]>0.1) score <- 2
    else if(temp[i,"AUC.Antagonist"]>0.001) score <- 1
    if(temp[i,"specificity_Z"]>0.8) score <- score+1
    score <- score+temp[i,"antagonist_confirmation_score"]
    temp[i,"antagonist_confidence_score"] <- score
  }
  SUPERMATRIX_3 <<- temp
  file <- "../output/superMatrix_3.xlsx"
  sty <- createStyle(halign="center",valign="center",textRotation=90,textDecoration = "bold")
  write.xlsx(SUPERMATRIX_3,file,firstActiveRow=2,firstActiveCol=5,headerStyle=sty)
  
}