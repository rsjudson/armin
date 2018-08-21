source("general_utils.R")
#--------------------------------------------------------------------------------------
#'
#' Prepare the cytotoxicity files
#'
#--------------------------------------------------------------------------------------
prepCytotox <- function() {
  printCurrentFunction()
   file <- "../input/AR_pathway_chemicals.xlsx"
  chems <- read.xlsx(file)
  code.list <- chems[,"code"]
  file <- "../input/cytotox_ranges_new.xlsx"
  cyto <- read.xlsx(file)
  rownames(cyto) <- cyto[,"code"]
  cyto <- cyto[code.list,]
  ratio <- cyto[,"nhit"]/(cyto[,"ntested"]+1)
  x <- cyto[,"cytotox_median_log"] 
  x[ratio<0.1] <- 2
  x[cyto[,"cytotox_median_log"]==3] <- 3
  cyto[,"cytotox_median_log"] <- x
  file <- "../input/AR_pathway_cytotox.xlsx"
  write.xlsx(cyto,file)
}