source("general_utils.R")
#--------------------------------------------------------------------------------------
#'
#' Prepare the input files
#'
#--------------------------------------------------------------------------------------
prepFlags <- function(date.string="180614") {
  printCurrentFunction()
  file <- "../input/assay_list.txt"
  assay.list <- readLines(file,n=-1)
  file <- "../input/AR_pathway_chemicals.xlsx"
  chems <- read.xlsx(file)
  code.list <- chems[,"code"]
  
  file <- paste0("../input/toxcast_matrix/toxcast_flags_",date.string,".csv")
  print(file)
  mat <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  mat <- mat[is.element(mat[,"code"],code.list),]
  mat <- mat[is.element(mat[,"aenm"],assay.list),]
  
  file <- "../input/AR_pathway_flags.xlsx"
  write.xlsx(mat,file)
}