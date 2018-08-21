#--------------------------------------------------------------------------------------
#'
#' load data required for the calculations
#'
#--------------------------------------------------------------------------------------
loadData <- function() {
  printCurrentFunction()

  file <- "../input/AR_pathway_chemicals.xlsx"
  chems <- read.xlsx(file)
  rownames(chems) <- chems[,"code"]
  CHEMS <<- chems
  CODE.LIST <<- chems[,"code"]
  NCHEM <<- length(CODE.LIST)
  cat("loaded Chemical Information\n")
  cat("Unique chemicals: ",NCHEM,"\n")
  
  file <- "../input/AR_pathway_logac50.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp2 <- temp[,5:dim(temp)[2]]
  MAT.logac50 <<- temp2
  ASSAY.LIST <<- names(MAT.logac50)
  
  file <- "../input/AR_pathway_logac50_cytotox_filtered.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp2 <- temp[,5:dim(temp)[2]]
  MAT.logac50.cytofiltered <<- temp2
 
  file <-"../input/AR_pathway_z.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.z <<- temp
  
  file <-"../input/AR_pathway_cytotox.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  CYTOTOX <<- temp
  
  file <- "../input/AR_pathway_t.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  MAT.t <<- temp
  
  file <- "../input/AR_pathway_w.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  MAT.w <<- temp
  
  file <- "../input/AR_pathway_emax.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.emax <<- temp
  
  file <- "../input/AR_pathway_max_conc.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.maxconc <<- 10**temp
  
  file <- "../input/AR_pathway_logac10.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.logac10 <<- temp
  
  file <- "../input/AR_pathway_logacb.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.logacb <<- temp
  
  file <- "../input/AR_pathway_logacc.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.logacc <<- temp
  
  file <- "../input/AR_pathway_hitcall.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  temp <- temp[,5:dim(temp)[2]]
  MAT.hitcall <<- temp

  file <- "../input/AR_pathway_refchems.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp[,"code"]
  REFCHEMS <<- temp
  
  file <- "../input/tox21_qc.txt"
  temp <- read.table(file,header=T,sep="\t",stringsAsFactors=F,quote="")
  TOX21.QC <<- temp
  
  file <- "../input/AR_pathway_flags.xlsx"
  temp <- read.xlsx(file)
  names(temp)[7] <- "code"
  names(temp)[5] <- "name"
  names(temp)[9] <- "assay"
  CAUTION.FLAGS <<- temp
  
  file <- "../input/EDSP_universe_input_2014_04_17.txt"
  EDSP.UNIVERSE <<- read.table(file,sep="\t",header=T,stringsAsFactors=F,comment.char="",quote="\"")
}
