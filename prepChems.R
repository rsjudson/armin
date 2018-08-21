source("general_utils.R")
library(openxlsx)
#--------------------------------------------------------------------------------------
#'
#' Prepare the chemical list
#'
#--------------------------------------------------------------------------------------
prepChems <- function(date.string="180614") {
  printCurrentFunction()
  file <- "../input/assay_list.txt"
  assay.list <- readLines(file,n=-1)
  file <- paste0("../input/toxcast_matrix/ac50_Matrix_",date.string,".csv")
  print(file)
  mat <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  rownames(mat) <- mat[,1]
  temp <- mat[,assay.list]
  temp[!is.na(temp)] <- 1
  temp[is.na(temp)] <- 0
  rs <- rowSums(temp)
  hist(rs)
  temp <- temp[rs==length(assay.list),]
  cat("rows left: ",dim(temp)[1],"\n")
  code.list <- rownames(temp)
  file <- "../input/toxcast_matrix/toxcast_chemicals_2018-08-21.xlsx"
  chems <- read.xlsx(file)
  rownames(chems) <- chems[,"code"]
  chems <- chems[code.list,]
  if(!exists("DSSTOX")){ 
    file <- "../input/DSStox/DSSTOX.RData"
    load(file=file)
    DSSTOX <<- DSSTOX
  }
  chems <- cbind(chems[,1],chems)
  names(chems)[1] <- "dsstox_substance_id"
  casrn.list <- chems[,"casrn"]
  dtx.list <- DSSTOX[is.element(DSSTOX[,"casrn"],casrn.list),"dsstox_substance_id"]
  chems[,1] <- dtx.list
  file <- "../input/AR_pathway_chemicals.xlsx"
  write.xlsx(chems,file)
  
}