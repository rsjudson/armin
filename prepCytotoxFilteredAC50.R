source("general_utils.R")
#--------------------------------------------------------------------------------------
#'
#' Prepare the cytotoxicity files
#'
#--------------------------------------------------------------------------------------
prepCytotoxFilteredAC50 <- function() {
  printCurrentFunction()
  file <- "../input/AR_pathway_cytotox.xlsx"
  cyto <- read.xlsx(file)
  file <- "../input/AR_pathway_logac50.xlsx"
  mat <- read.xlsx(file)
  cutoff <- log10(cyto[,"cytotox_lower_bound_um"]) # gives uM value
  file <- "../input/assay_list.txt"
  assay.list <- readLines(file,n=-1)
  for(assay in assay.list) {
    x <- mat[,assay]
    x[x>cutoff] <- 0
    mat[,assay] <- x
  }
  file <- "../input/AR_pathway_logac50_cytotox_filtered.xlsx"
  write.xlsx(mat,file)
}