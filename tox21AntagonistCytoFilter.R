source("general_utils.R")
#--------------------------------------------------------------------------------------
#'
#' Prepare the cytotoxicity files
#'
#--------------------------------------------------------------------------------------
tox21AntagonistCytoFilter <- function(date.string="180614") {
  printCurrentFunction()
  file <- "../input/AR_pathway_logac50.xlsx"
  ac50 <- read.xlsx(file)
  code.list <- ac50[,"code"]
  file <- paste0("../input/toxcast_matrix/modl_ga_Matrix_",date.string,".csv")
  toxcast <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  rownames(toxcast) <- toxcast[,1]
  toxcast <- toxcast[code.list,]
  
  file <- paste0("../input/toxcast_matrix/hitc_Matrix_",date.string,".csv")
  hitc <- read.table(file,stringsAsFactors=F,header=T,sep=",")
  rownames(hitc) <- hitc[,1]
  hitc <- hitc[code.list,]
  
  
  assay1 <- "TOX21_AR_BLA_Antagonist_ratio"  
  assay2 <- "TOX21_AR_BLA_Antagonist_viability" 
  x <- ac50[,assay1]
  y <- toxcast[,assay2]
  y[is.na(y)] <- 6
  z <- hitc[,assay2]
  z[is.na(z)] <- 0
  y[z==0] <- 6
  x[x>y] <- 6
  ac50[,assay1] <- x
  assay1 <- "TOX21_AR_LUC_MDAKB2_Antagonist2"
  assay2 <- "TOX21_AR_LUC_MDAKB2_Antagonist2_viability"
  x <- ac50[,assay1]
  y <- toxcast[,assay2]
  y[is.na(y)] <- 6
  z <- hitc[,assay2]
  z[is.na(z)] <- 0
  y[z==0] <- 6
  x[x>y] <- 6
  ac50[,assay1] <- x
  file <- "../input/AR_pathway_logac50.xlsx"
  write.xlsx(ac50,file)
}