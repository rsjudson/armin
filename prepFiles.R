source("general_utils.R")
library(openxlsx)
#--------------------------------------------------------------------------------------
#'
#' Prepare the input files
#' 
#' All concentrations are in log10(uM)
#'
#--------------------------------------------------------------------------------------
prepFiles <- function(date.string="180821",fix.nvs=F) {
  printCurrentFunction()
  file <- "../input/assay_list.txt"
  assay.list <- readLines(file,n=-1)
  file <- "../input/AR_pathway_chemicals.xlsx"
  chems <- read.xlsx(file)
  code.list <- chems[,"code"]
  
  file <- paste0("../input/toxcast_matrix/hitc_Matrix_",date.string,".csv")
  print(file)
  mat <- read.table(file,stringsAsFactors=F,header=T,sep=",")

  rownames(mat) <- mat[,1]
  mat <- mat[code.list,assay.list]
  mat[is.na(mat)] <- 0
  mat.hitcall <- mat
  mat <- cbind(chems,mat)
  file <- paste0("../input/AR_pathway_hitcall.xlsx")
  write.xlsx(mat,file)
  
  class.list <-   c("logac50","logac10",  "logacc",  "logacb",  "emax",    "t",      "w",      "z",     "max_conc")
  prefix.list <-  c("modl_ga","modl_ac10","modl_acc","modl_acb","resp_max","modl_tp","modl_gw","zscore","logc_max")
  default.list <- c(6,        6,          6,         6,          0,         0,        1,        0,        2)
  for(i in 1:length(class.list)){
    class <- class.list[i]
    prefix <- prefix.list[i]
    default <- default.list[i]
    
    file <- paste0("../input/toxcast_matrix/",prefix,"_Matrix_",date.string,".csv")
    print(file)
    mat <- read.table(file,stringsAsFactors=F,header=T,sep=",")
    rownames(mat) <- mat[,1]
    mat <- mat[code.list,assay.list]
    mat[is.na(mat)] <- default
    
    # fix the units for the NVS data
    if(fix.nvs) {
      if(is.element(class,c("logac50","logac10",  "logacc",  "logacb", "max_conc"))) {
        assay.list.nvs <- c("NVS_NR_hAR","NVS_NR_rAR","NVS_NR_cAR")
        for(assay in assay.list.nvs) {
          x <- mat[,assay]
          if(class=="max_conc") {
            mat[,assay] <- default
          }
          else {
            x <- x+6
            x[x>=12] <- 6
            x[x>6] <- 6
            mat[,assay] <- x
          }
          #cat(assay,range(x),"\n")
        }
      }
    }
    #browser()
    if(class=="t") {
      mat[,"ATG_AR_TRANS_up"] <- mat[,"ATG_AR_TRANS_up"]*25
      mat[,"ACEA_AR_antagonist_80hr"] <- mat[,"ACEA_AR_antagonist_80hr"]*30
    }
    mat[mat.hitcall==0] <- default
    mat <- cbind(chems,mat)
    file <- paste0("../input/AR_pathway_",class,".xlsx")
    write.xlsx(mat,file)
  }
}