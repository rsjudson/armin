#--------------------------------------------------------------------------------------
#
# add the antagonist shift
#
#--------------------------------------------------------------------------------------
addAntagonistShift <- function(smatrix=SUPERMATRIX_1) {
  printCurrentFunction()
  temp <- smatrix
  name.list <- c("ga_min_left","ga_med_left","ga_max_left","hit_pct_left","ga_min_right","ga_med_right","ga_max_right","hit_pct_right")
  temp2 <- as.data.frame(matrix(nrow=dim(temp)[1],ncol=length(name.list)))
  names(temp2) <- name.list
  temp <- cbind(temp,temp2)

  file <- "../input/Tox21ARAntagonist_mc7.xlsx"
  mc7 <- read.xlsx(file)
  
  temp[,"antagonist_confirmation_score"] <- 0

  nchem <- dim(temp)[1]
  for(i in 1:nchem) {
    code <- temp[i,"code"]
    casrn <- temp[i,"casrn"]
    mtemp <- mc7[is.element(mc7[,"casn"],casrn),]
    score <- 0
    if(dim(mtemp)[1]>0) {
      temp_left <- mtemp[is.element(mtemp[,"assay_component_endpoint_name"],"TOX21_AR_LUC_MDAKB2_Antagonist2"),]
      if(dim(temp_left)[1]>0) {
        if(sum(temp_left[,"total_hitc"])==0) temp_left <- temp_left[1,]
        else {
          temp_left <- temp_left[order(temp_left[,"hit_pct"],decreasing=T),]
          temp_left <- temp_left[1,]
        }
      }
      
      temp_right <- mtemp[is.element(mtemp[,"assay_component_endpoint_name"],"TOX21_AR_LUC_MDAKB2_Antagonist"),]
      if(dim(temp_right)[1]>0) {
        if(sum(temp_right[,"total_hitc"])==0) temp_right <- temp_right[1,]
        else {
          temp_right <- temp_right[order(temp_right[,"hit_pct"],decreasing=T),]
          temp_right <- temp_right[1,]
        }
      }
      
      if(temp_left[1,"hit_pct"]>0.5 && temp_right[1,"hit_pct"]>0.5) {
        if(temp_left[1,"modl_ga_med"]<temp_right[1,"modl_ga_med"]) {
          if(temp_left[1,"modl_ga_max"]<temp_right[1,"modl_ga_min"]) score <- 3
          else score <- 1
        }
        else score <- -1
      }
      else {
        if(temp_left[1,"hit_pct"]>0.5 && temp_right[1,"hit_pct"]<0.5) score <- 2
        else if(temp_left[1,"hit_pct"]<0.5 && temp_right[1,"hit_pct"]>0.5) score <- -1
      }
    }
    temp[i,"ga_min_left"] <- temp_left[1,"modl_ga_min"]
    temp[i,"ga_med_left"] <- temp_left[1,"modl_ga_med"]
    temp[i,"ga_max_left"] <- temp_left[1,"modl_ga_max"]
    temp[i,"hit_pctc_left"] <- temp_left[1,"hit_pct"]
    temp[i,"ga_min_right"] <- temp_right[1,"modl_ga_min"]
    temp[i,"ga_med_right"] <- temp_right[1,"modl_ga_med"]
    temp[i,"ga_max_right"] <- temp_right[1,"modl_ga_max"]
    temp[i,"hit_pct_right"] <- temp_right[1,"hit_pct"]
    temp[i,"antagonist_confirmation_score"] <- score
    #if(casrn=="84371-65-3") browser()
  }
  SUPERMATRIX_2 <<- temp
  file <- "../output/superMatrix_2.xlsx"
  sty <- createStyle(halign="center",valign="center",textRotation=90,textDecoration = "bold")
  write.xlsx(SUPERMATRIX_2,file,firstActiveRow=2,firstActiveCol=4,headerStyle=sty)
}