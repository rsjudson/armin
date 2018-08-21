#--------------------------------------------------------------------------------------
#
# Calculate a composite approximate AC50 for a single chemical
#
#===============================================================================
# Pathway-specific - start
#===============================================================================
#------------------------------------------------------------------------------------
pseudoAC <- function(code,variable="AC50",AUC.R1,AUC.R2,debug=F) {
  cutoff <- 0.05
  
  agonist.assays <- c("NVS_NR_hAR","NVS_NR_cAR","NVS_NR_rAR","OT_AR_ARSRC1_0480","OT_AR_ARSRC1_0960","ATG_AR_TRANS_up","OT_AR_ARELUC_AG_1440","TOX21_AR_BLA_Agonist_ratio","TOX21_AR_LUC_MDAKB2_Agonist","ACEA_AR_agonist_80hr")
  antagonist.assays <- c("NVS_NR_hAR","NVS_NR_cAR","NVS_NR_rAR","OT_AR_ARSRC1_0480","OT_AR_ARSRC1_0960","TOX21_AR_BLA_Antagonist_ratio","TOX21_AR_LUC_MDAKB2_Antagonist2","ACEA_AR_antagonist_80hr")
  
  median.agonist <- 6
  median.antagonist <- 6
  min.agonist <- 6
  min.antagonist <- 6
  if(AUC.R1>cutoff) {
    if(variable=="AC50") temp <- MAT.logac50[code,agonist.assays]
    if(variable=="AC10") temp <- MAT.logac10[code,agonist.assays]
    if(variable=="ACC") temp <- MAT.logacc[code,agonist.assays]
    if(variable=="ACB") temp <- MAT.logacb[code,agonist.assays]
    temp <- temp[!is.na(temp)]
    temp <- temp[!is.nan(temp)]
    temp <- temp[temp<6]
    if(length(temp)>0) {
      tmedian <- median(temp)
      if(length(temp)>2) tmad <- mad(temp)
      else tmad <- 0.5
      median.agonist <- tmedian
      min.agonist <- tmedian-3*tmad
    }
  }
  if(AUC.R2>cutoff) {
    if(variable=="AC50") temp <- MAT.logac50[code,antagonist.assays]
    if(variable=="AC10") temp <- MAT.logac10[code,antagonist.assays]
    if(variable=="ACC") temp <- MAT.logacc[code,antagonist.assays]
    if(variable=="ACB") temp <- MAT.logacb[code,antagonist.assays]
    temp <- temp[!is.na(temp)]
    temp <- temp[!is.nan(temp)]
    temp <- temp[temp<6]
    if(length(temp)>0) {
      tmedian <- median(temp)
      if(length(temp)>2) tmad <- mad(temp)
      else tmad <- 0.5
      median.antagonist <- tmedian
      min.antagonist <- tmedian-3*tmad
    }
  }
  value <- list(median.value=min(median.agonist,median.antagonist),
                min.value=min(min.agonist,min.antagonist))
  if(debug){ 
    cat("AUC.R1:",AUC.R1,"\n")
    cat("AUC.R2:",AUC.R2,"\n")
    print(value)
    browser()
  }
  return(value)
}
