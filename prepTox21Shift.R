source("general_utils.R")
#--------------------------------------------------------------------------------------
#'
#' Prepare the input files
#' 
#' All concentrations are in log10(uM)
#'
#--------------------------------------------------------------------------------------
prepTox21Shift <- function(date.string="180409") {
  printCurrentFunction()
  code.list <- CODE.LIST
  file <- paste0("../input/toxcast_matrix/modl_ga_Matrix_",date.string,".csv")
  toxcast <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  rownames(toxcast) <- toxcast[,1]
  toxcast <- toxcast[code.list,]
  
  file <- paste0("../input/toxcast_matrix/hitc_Matrix_",date.string,".csv")
  hitmat <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  rownames(hitmat) <- hitmat[,1]
  hitmat <- hitmat[code.list,]
 
  assay1 <- "TOX21_AR_LUC_MDAKB2_Antagonist2"
  assay2 <- "TOX21_AR_LUC_MDAKB2_Antagonist2_viability"
  assy <- toxcast[,assay1]
  assy[is.na(assy)] <- 6
  hitc <- hitmat[,assay1]
  hitc[is.na(hitc)] <- 0
  assy[hitc==0] <- 6
  
  cyto <- toxcast[,assay2]
  cyto[is.na(cyto)] <- 6
  hitc <- hitmat[,assay2]
  hitc[is.na(hitc)] <- 0
  cyto[hitc==0] <- 6
  
  assy[assy>cyto] <- 6
  antagonist2 <- assy
  
  assay1 <- "TOX21_AR_LUC_MDAKB2_Antagonist"
  assay2 <- "TOX21_AR_LUC_MDAKB2_Antagonist_viability"
  assy <- toxcast[,assay1]
  assy[is.na(assy)] <- 6
  hitc <- hitmat[,assay1]
  hitc[is.na(hitc)] <- 0
  assy[hitc==0] <- 6
  
  cyto <- toxcast[,assay2]
  cyto[is.na(cyto)] <- 6
  hitc <- hitmat[,assay2]
  hitc[is.na(hitc)] <- 0
  cyto[hitc==0] <- 6
  
  assy[assy>cyto] <- 6
  antagonist1 <- assy
  
  delta <- antagonist2 - antagonist1
  name.list <- c(names(CHEMS),"TOX21_AR_LUC_MDAKB2_Antagonist","TOX21_AR_LUC_MDAKB2_Antagonist2","delta")
  mat <- as.data.frame(matrix(nrow=length(CODE.LIST),ncol=length(name.list)))
  names(mat) <- name.list
  mat[,1:4] <- CHEMS
  mat[,"TOX21_AR_LUC_MDAKB2_Antagonist"] <- antagonist1
  mat[,"TOX21_AR_LUC_MDAKB2_Antagonist2"] <- antagonist2
  mat[,"delta"] <- delta
  
  file <- paste0("../input/AR_pathway_tox21_logac50_shift.xlsx")
  write.xlsx(mat,file)
}