#--------------------------------------------------------------------------------------
#
# ER_assay_select_v02.R - select a minimal set of assays for the ER model
#
# August 2016
# Richard Judson
#
# US EPA
# Questions, comments to: judson.richard@epa.gov, 919-541-3085
#
# File structure: 
#
# bootstrap: input file from the bootstrap sampling appraoch of Watt et al. 
# manuscript: copies of the manuscript and supplemental files
# output: Excel files produced by the code
# plots: pdf plots produced by the code
# refchems: input files concerning reference chemicals: Table 2
# R: code
#
# Steps to run all code:
#
# 1. prep.data.2() - prepares all of data frames in memory
# 2. export.AUC() - exports Supplemental table S1.med, .max, .min
# 3. cutoff.scan(T,2) - runs the scan to determine the best cutoff, and provides evidence: ends up being 0.1
# 4. assay.scan.2(T,"allchems",0.1,2) - run the scan across all combinations of 1 to 16 assays, Figure 1
# 5. summarize.sample.stats() - output the stats for all subset models
# 6. plot.stats.1(T) - Figure 2
# 7. index.sample.stats.1(T,T) more summary of the statistics, Figure 3, Table 3
# 8. export.all.data(T)  - yields a spreadsheet with all models above the specified threshold (supplemental file S2), plus Table 4
# 9. count.misclassified.chems(T) - generates informaiton for Figure 4
# 10. misclass.hm(T) - Figure 4
#
#--------------------------------------------------------------------------------------
library(grDevices)
library(RColorBrewer)
library(stringr)
library(openxlsx)
source("utils_20170321.R")

#--------------------------------------------------------------------------------------
#
# prepare the data set
#
#--------------------------------------------------------------------------------------
prep.data.2 <- function() {
  print.current.function ()
  #  assay.list <- c("NVS_NR_bER","NVS_NR_hER","NVS_NR_mERa","OT_ER_ERaERa_0480","OT_ER_ERaERa_1440","OT_ER_ERaERb_0480","OT_ER_ERaERb_1440","OT_ER_ERbERb_0480","OT_ER_ERbERb_1440","OT_ERa_EREGFP_0120","OT_ERa_EREGFP_0480","ATG_ERa_TRANS_up","ATG_ERE_CIS_up","Tox21_ERa_BLA_Agonist_ratio","Tox21_ERa_LUC_BG1_Agonist","ACEA_T47D_80hr_Positive","Tox21_ERa_BLA_Antagonist_ratio","Tox21_ERa_LUC_BG1_Antagonist")
  assay.list <- c("NVS_NR_bER","NVS_NR_hER","NVS_NR_mERa","OT_ER_ERaERa_0480","OT_ER_ERaERa_1440","OT_ER_ERaERb_0480","OT_ER_ERaERb_1440","OT_ER_ERbERb_0480","OT_ER_ERbERb_1440","OT_ERa_EREGFP_0120","OT_ERa_EREGFP_0480","ATG_ERa_TRANS_up","ATG_ERE_CIS_up","TOX21_ERa_BLA_Agonist_ratio","TOX21_ERa_LUC_BG1_Agonist","ACEA_T47D_80hr_Positive")
  assay.class <- c("R3","R3","R3","R4","R4","R4","R4","R4","R4","R5","R5","R6","R6","R7","R7","R8")

  file <- "../bootstrap/er_auc_results_for_plot.rds"
  full.model <- readRDS(file)
  full.model <- full.model[is.element(full.model[,"auc_type"],"AUC.Agonist"),]
  temp <- full.model[,c("code","chnm","auc_value","minbar","maxbar")]
  names(temp) <- c("CODE","Name","AUC.agonist","minbar","maxbar")
  rownames(temp) <- temp[,"CODE"]
  FULL.MODEL <<- temp
  
  file <- "../bootstrap/er_model_assay_subset_predictivity.rds"
  bootstrap.data <- readRDS(file)
  BOOTSTRAP.DATA <<- bootstrap.data
  ASSAY.LIST <<- assay.list
  ASSAY.CLASS <<- assay.class
  NASSAY <<- length(ASSAY.LIST)
  NCHEM <<- dim(FULL.MODEL)[1]
  fix.assay.auc()
  prep.refchems()
  summarize.boot()
  normalize.auc()
}
#--------------------------------------------------------------------------------------
#
# calculate the corrected assay-wise AUC values
#
#--------------------------------------------------------------------------------------
fix.assay.auc <- function() {
  print.current.function ()
  boot <- BOOTSTRAP.DATA
  boot <- boot[is.element(boot[,"aenm"],ASSAY.LIST),]
  temp <- boot[,c("code","aenm","auc")]
  top <- boot[,"modl_tp"]
  ga <- boot[,"modl_ga"]
  top[is.na(top)] <- 0
  ga[is.na(ga)] <- 0
  top[top>100] <- 100
  top <- top/100
  ga <- 6-ga
  auc <- top*ga
  temp[,"auc"] <- auc
  names(temp) <- c("CODE","assay","AUC")
  BOOT.AUC <<- temp
}
#--------------------------------------------------------------------------------------
#
# Prep the reference chemical information
#
#--------------------------------------------------------------------------------------
prep.refchems <- function() {
  print.current.function ()
  file <- "../refchems/NICEATM_refchems.xlsx"
  refchems <- read.xlsx(file)
  refchems <- refchems[refchems[,"useme"]==1,]
  code.list <- unique(sort(refchems[,"CODE"]))
  name.list <- c("CODE","CASRN","Name","InVitro","InVivo")
  nchem <- length(code.list)
  mat <- as.data.frame(matrix(nrow=nchem,ncol=length(name.list)))
  names(mat) <- name.list
  mat[,"CODE"] <- code.list
  rownames(mat) <- code.list
  for(i in 1:nchem) {
    code <- code.list[i]
    temp <- refchems[is.element(refchems[,"CODE"],code),]
    for(j in 1:dim(temp)[1]) {
      mat[i,"Name"] <- temp[j,"Name"]
      mat[i,"CASRN"] <- temp[j,"CASRN"]
      aclass <- temp[j,"class"]
      #print(temp[j,])
      if(aclass=="in vivo") {
        mat[i,"InVivo"] <- temp[j,"in.vivo.call"]
      }
      else {
        mat[i,"InVitro"] <- temp[j,"in.vitro.strength"]
      }
    }
  }
  file <- "../refchems/refchem_table.xlsx"
  write.xlsx(mat,file)
  for(i in 1:nchem) {
    if(!is.na(mat[i,"InVivo"])) {
      if(mat[i,"InVivo"]=="Active") mat[i,"InVivo"] <- 1
      if(mat[i,"InVivo"]=="Inactive") mat[i,"InVivo"] <- 0
    }
    if(!is.na(mat[i,"InVitro"])) {
      if(mat[i,"InVitro"]=="Inactive") mat[i,"InVitro"] <- 0
      else mat[i,"InVitro"] <- 1
    }
  }
  REFCHEM.DATA <<- mat
  REFCHEM.CODES <<- code.list
}
#--------------------------------------------------------------------------------------
#
# calculate the median CI values for all chemicals, all assays
#
#--------------------------------------------------------------------------------------
summarize.boot <- function() {
  print.current.function ()
  temp <- as.data.frame(matrix(nrow=NCHEM,ncol=NASSAY))
  code.list <- rownames(FULL.MODEL)
  rownames(temp) <- code.list
  names(temp) <- ASSAY.LIST
  temp[] <- 0
  mat.med <- temp
  mat.min <- temp
  mat.max <- temp
  for(i in 1:NCHEM) {
    code <- code.list[i]
    cname <- FULL.MODEL[code,"Name"]
    cat(code,":",cname,"\n")
    temp <- BOOT.AUC[is.element(BOOT.AUC[,"CODE"],code),]
    assay.list <- unique(temp[,"assay"])
    nassay <- length(assay.list)
    for(j in 1:nassay) {
      assay <- assay.list[j]
      x <- temp[is.element(temp[,"assay"],assay),"AUC"]
      if(max(x)>0) {
        mat.med[code,assay] <- median(x)
        q <- quantile(x,probs=seq(0,1,0.025))
        mat.min[code,assay] <- q[2]
        mat.max[code,assay] <- q[40]
      }
    }
  }
  MAT.MED <<- mat.med
  MAT.MAX <<- mat.max
  MAT.MIN <<- mat.min
  file <- "../bootstrap/mat.med.rds"
  saveRDS(mat.med,file)
  file <- "../bootstrap/mat.min.rds"
  saveRDS(mat.min,file)
  file <- "../bootstrap/mat.max.rds"
  saveRDS(mat.max,file)
}
#--------------------------------------------------------------------------------------
#
# normalize the AUC matices
#
#--------------------------------------------------------------------------------------
normalize.auc <- function() {
  print.current.function ()
  file <- "../bootstrap/mat.med.rds"
  mat.med <- readRDS(file)
  file <- "../bootstrap/mat.min.rds"
  mat.min <- readRDS(file)
  file <- "../bootstrap/mat.max.rds"
  mat.max <- readRDS(file)
  
  maxval <- colMax(mat.med)
  nassay <- dim(mat.med)[2]
  for(i in 1:nassay) {
    mat.med[,i] <- mat.med[,i]/maxval[i]
    mat.min[,i] <- mat.min[,i]/maxval[i]
    mat.max[,i] <- mat.max[,i]/maxval[i]
  }
  MAT.MED <<- mat.med
  MAT.MAX <<- mat.max
  MAT.MIN <<- mat.min
  file <- "../bootstrap/mat.med.rds"
  saveRDS(mat.med,file)
  file <- "../bootstrap/mat.min.rds"
  saveRDS(mat.min,file)
  file <- "../bootstrap/mat.max.rds"
  saveRDS(mat.max,file)
}
#--------------------------------------------------------------------------------------
#
# Export the normalized median AUC matrices
#
#--------------------------------------------------------------------------------------
export.AUC <- function() {
  print.current.function ()
  x <- cbind(FULL.MODEL,MAT.MED)
  file <- "../output/S1a_auc_med.xlsx"
  write.xlsx(x,file)
  x <- cbind(FULL.MODEL,MAT.MIN)
  file <- "../output/S1b_auc_min.xlsx"
  write.xlsx(x,file)
  x <- cbind(FULL.MODEL,MAT.MAX)
  file <- "../output/S1c_auc_max.xlsx"
  write.xlsx(x,file)
}
#--------------------------------------------------------------------------------------
#
# cutoff.scan
#
#--------------------------------------------------------------------------------------
cutoff.scan <- function(to.file=F,minhit=2) {
  print.current.function ()
  if(to.file) {
    fname <- paste("../plots/Fig 1 cutoff_scan_",minhit,".pdf",sep="")
    pdf(file=fname,width=5.5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,4))

  code.list <- FULL.MODEL[,"CODE"]
  code.list.refchems <- REFCHEM.CODES
  nchem.ref <- length(code.list.refchems)
  code.list.nonref <- code.list[!is.element(code.list,code.list.refchems)]
  nchem.nonref <- length(code.list.nonref)
  
  auc.model.refchems <- FULL.MODEL[code.list.refchems,]
  auc.model.nonref   <- FULL.MODEL[code.list.nonref,]
  auc.model <- rbind(auc.model.refchems,auc.model.nonref)
  nchem <- dim(auc.model)[1]
  
  name.list <- c("cutoff","RSE.allchems","R2.allchems","RSE.refchems","R2.refchems","sens.allchem","spec.allchem","ba.allchem","sens.refchem","spec.refchem","ba.refchem")
  res <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  res.all <- NULL
  names(res) <- name.list

  code.list <- c(code.list.refchems,code.list.nonref)
  mat.med <- MAT.MED[code.list,]
  mat.min <- MAT.MIN[code.list,]
  mat.max <- MAT.MAX[code.list,]
  bamax <- 0
  amask <- vector(length=NASSAY,mode="integer")
  amask[] <- 1
  cutoff.list <- c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.04,0.03,0.05,0.06,0.07,0.08,0.09,0.1)
  ncutoff <- length(cutoff.list)
  for(i in 1:ncutoff) {
    cutoff <- cutoff.list[i]
    cat("==========================================================\n")
    cat(cutoff,"\n")
    cat("==========================================================\n")
    temp <- mat.med[,amask==1]
    temp[temp>0] <- 1
    rs <- rowSums(temp)
    hit.filter <- rs
    hit.filter[rs<minhit] <- 0
    hit.filter[hit.filter>0] <- 1
    if(sum(amask)<minhit) hit.filter[] <- 1
    x <- as.matrix(mat.med[,amask==1] * hit.filter)
    y <- auc.model[,"AUC.agonist"] * hit.filter
    
    mask <- y
    mask[mask>cutoff] <- 1
    x <- x[mask==1,]
    y <- y[mask==1]
    model <- lm(y ~ x + 0)
    print(summary(model))
    coefs <- model$coefficients
    rse.allchems <- summary(model)$sigma
    r2.allchems <- summary(model)$adj.r.squared
    res[1,"cutoff"] <- cutoff
    res[1,"RSE.allchems"] <- rse.allchems
    res[1,"R2.allchems"] <- r2.allchems
    
    xref <- as.matrix(mat.med[1:nchem.ref,amask==1] * hit.filter[1:nchem.ref])
    yref <- auc.model[1:nchem.ref,"AUC.agonist"] * hit.filter[1:nchem.ref]
    
    model.ref <- lm(yref ~ xref+0)
    rse.refchems <- summary(model.ref)$sigma
    r2.refchems <- summary(model.ref)$adj.r.squared
    
    res[1,"RSE.refchems"] <- rse.refchems
    res[1,"R2.refchems"] <- r2.refchems
    
    xp <- as.matrix(mat.med[,amask==1]) %*% coefs
    xp.min <- as.matrix(mat.min[,amask==1]) %*% coefs
    xp.max <- as.matrix(mat.max[,amask==1]) %*% coefs
    xp[xp<0] <- 0
    xp.min[xp.min<0] <- 0
    xp.max[xp.max<0] <- 0
    
    yp <- auc.model[,"AUC.agonist"]
    xp2 <- xp
    yp2 <- yp
    xp2[xp2<cutoff] <- 0
    xp2[xp2>0] <- 1
    yp2[yp2<cutoff] <- 0
    yp2[yp2>0] <- 1
    a <- sum(xp2*yp2)
    b <- sum(xp2*(1-yp2))
    c <- sum((1-xp2)*yp2)
    d <- sum((1-xp2)*(1-yp2))
    txt <- TxT(a,b,c,d)
    res[1,"sens.allchem"] <- txt$sens
    res[1,"spec.allchem"] <- txt$spec
    res[1,"ba.allchem"] <- txt$ba
    ba.all <- txt$ba

    main <- paste("cutoff: ",cutoff)
    plot(yp~xp,xlim=c(0,1),ylim=c(0,1),xlab="Subset Model AUC Score",ylab="Full Model AUC Score",main=main)
    for(j in 1:nchem) {
      x <- xp[j]
      y <- yp[j]
      minval <- auc.model[j,"minbar"]
      maxval <- auc.model[j,"maxbar"]
      lines(c(x,x),c(minval,maxval))
      if(x>=minval && x<=maxval) points(x,y,bg="red",cex=1,pch=21)
      minval <- xp.min[j]
      maxval <- xp.max[j]
      lines(c(minval,maxval),c(y,y))
      if(y>=minval && y<=maxval) points(x,y,bg="red",cex=1,pch=21)
    }      
    lines(c(0,1),c(0,1))
    lines(c(0,1),c(cutoff,cutoff))
    lines(c(cutoff,cutoff),c(0,1))
    text(cutoff,0.99,paste("Assays:",sum(amask)),pos=4,cex=0.9)
    text(cutoff,0.92,paste("RMSE all chems:",format(rse.allchems,digits=2)),pos=4,cex=0.9)
    text(cutoff,0.85,paste("R2 all chems:",format(r2.allchems,digits=2)),pos=4,cex=0.9)
    text(cutoff,0.78,paste("RMSE ref chems:",format(rse.refchems,digits=2)),pos=4,cex=0.9)
    text(cutoff,0.71,paste("R2 ref chems:",format(r2.refchems,digits=2)),pos=4,cex=0.9)
    text(0.5,cutoff+0.21,"Sens,Spec,BA",pos=4,cex=1)
    text(0.5,cutoff+0.14,paste("All: [",format(txt$sens,digits=2),":",format(txt$spec,digits=2),":",format(txt$ba,digits=2),"]",sep=""),pos=4,cex=1)
    
    xp2 <- xp2[1:nchem.ref]
    yp2 <- yp2[1:nchem.ref]
    a <- sum(xp2*yp2)
    b <- sum(xp2*(1-yp2))
    c <- sum((1-xp2)*yp2)
    d <- sum((1-xp2)*(1-yp2))
    txt <- TxT(a,b,c,d)      
    res[1,"sens.refchem"] <- txt$sens
    res[1,"spec.refchem"] <- txt$spec
    res[1,"ba.refchem"] <- txt$ba
    text(0.5,cutoff+0.07,paste("Ref: [",format(txt$sens,digits=2),":",format(txt$spec,digits=2),":",format(txt$ba,digits=2),"]",sep=""),pos=4,cex=1)
    
    # plot without labels
    main <- paste("cutoff: ",cutoff)
    plot(yp~xp,xlim=c(0,1),ylim=c(0,1),xlab="Subset Model AUC Score",ylab="Full Model AUC Score",main=main)
    for(j in 1:nchem) {
      x <- xp[j]
      y <- yp[j]
      minval <- auc.model[j,"minbar"]
      maxval <- auc.model[j,"maxbar"]
      lines(c(x,x),c(minval,maxval))
      if(x>=minval && x<=maxval) points(x,y,bg="red",cex=1,pch=21)
      minval <- xp.min[j]
      maxval <- xp.max[j]
      lines(c(minval,maxval),c(y,y))
      if(y>=minval && y<=maxval) points(x,y,bg="red",cex=1,pch=21)
    }      
    lines(c(0,1),c(0,1))
    lines(c(0,1),c(cutoff,cutoff))
    lines(c(cutoff,cutoff),c(0,1))
    text(cutoff,0.99,paste("Assays:",sum(amask)),pos=4,cex=0.9)

    xp2 <- xp2[1:nchem.ref]
    yp2 <- yp2[1:nchem.ref]
    a <- sum(xp2*yp2)
    b <- sum(xp2*(1-yp2))
    c <- sum((1-xp2)*yp2)
    d <- sum((1-xp2)*(1-yp2))
    txt <- TxT(a,b,c,d)      
    res[1,"sens.refchem"] <- txt$sens
    res[1,"spec.refchem"] <- txt$spec
    res[1,"ba.refchem"] <- txt$ba

    res.all <- rbind(res.all,res)
    if(!to.file) browser()
    
  }
  file <- "../output/cutoff.scan.xlsx"
  write.xlsx(res.all,file)
  if(to.file) dev.off()
}
#--------------------------------------------------------------------------------------
#
# assay.scan - produces Figure 1 of the paper
#
#--------------------------------------------------------------------------------------
assay.scan.2 <- function(to.file=F,chemset="allchems",cutoff=0.1,minhit=2) {
  print.current.function ()
  if(to.file) {
    fname <- paste("../plots/assay_scan_2_",chemset,"_",cutoff,".pdf",sep="")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,4,4,4))
  
  code.list <- FULL.MODEL[,"CODE"]

  code.list.refchems <- REFCHEM.CODES
  nchem.ref <- length(code.list.refchems)
  code.list.nonref <- code.list[!is.element(code.list,code.list.refchems)]
  nchem.nonref <- length(code.list.nonref)

  auc.model.refchems <- FULL.MODEL[code.list.refchems,]
  auc.model.nonref   <- FULL.MODEL[code.list.nonref,]
  auc.model <- rbind(auc.model.refchems,auc.model.nonref)
  nchem <- dim(auc.model)[1]

  scan.grid <- expand.grid(seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1))
  names(scan.grid) <- ASSAY.LIST
  rs <- rowSums(scan.grid)
  index <- sort(rs,index.return=T)$ix
  scan.grid <- scan.grid[index,]
  scan.grid <- scan.grid[2:dim(scan.grid)[1],]
  nscan <- dim(scan.grid)[1]
  
  name.list <- c("Assays","RSE.allchems","R2.allchems","RSE.refchems","R2.refchems",
                 "sens.allchem","spec.allchem","ba.allchem",
                 "sens.refchem","spec.refchem","ba.refchem",
                 "sens.refchem.invitro","spec.refchem.invitro","ba.refchem.invitro",
                 "sens.refchem.invivo","spec.refchem.invivo","ba.refchem.invivo",
                 "ba.min"
  )
  res <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(res) <- name.list
  res.all <- NULL
  res.amask <- as.data.frame(matrix(nrow=1,ncol=length(ASSAY.LIST)))
  names(res.amask) <- ASSAY.LIST
  res.amask.all <- NULL
  res.coefs <- as.data.frame(matrix(nrow=1,ncol=length(ASSAY.LIST)))
  for(i in 1:NASSAY) names(res.coefs)[i] <- paste(ASSAY.LIST[i],"_weight",sep="")
  res.coefs.all <- NULL
  
  ndo <- nscan
  time0 <- Sys.time()
  counter <- 0

  code.list <- c(code.list.refchems,code.list.nonref)
  mat.med <- MAT.MED[code.list,]
  mat.min <- MAT.MIN[code.list,]
  mat.max <- MAT.MAX[code.list,]
  bamax <- 0
  
  for(i in 1:ndo) {
    amask <- as.numeric(scan.grid[i,])
    temp <- mat.med[,amask==1]
    temp[temp>0] <- 1
    if(sum(amask)>1) rs <- rowSums(temp)
    else rs <- temp
    hit.filter <- rs
    hit.filter[rs<minhit] <- 0
    hit.filter[hit.filter>0] <- 1
    if(sum(amask)<minhit) hit.filter[] <- 1
    x <- as.matrix(mat.med[,amask==1] * hit.filter)
    y <- auc.model[,"AUC.agonist"] * hit.filter

    mask <- y
    mask[mask>cutoff] <- 1
    x <- x[mask==1,]
    y <- y[mask==1]
    model <- lm(y ~ x + 0)
    coefs <- model$coefficients
    rse.allchems <- summary(model)$sigma
    r2.allchems <- summary(model)$adj.r.squared
    res.amask[1,] <- amask
    res.coefs[1,] <- 0
    counter <- 0
    for(j in 1:NASSAY) {
      if(amask[j]==1) {
        counter <- counter+1
        res.coefs[1,j] <- coefs[counter]
      }
    }
    res[1,"Assays"] <- sum(amask)
    res[1,"RSE.allchems"] <- rse.allchems
    res[1,"R2.allchems"] <- r2.allchems
    
    xref <- as.matrix(mat.med[1:nchem.ref,amask==1] * hit.filter[1:nchem.ref])
    yref <- auc.model[1:nchem.ref,"AUC.agonist"] * hit.filter[1:nchem.ref]

    model.ref <- lm(yref ~ xref+0)
    rse.refchems <- summary(model.ref)$sigma
    r2.refchems <- summary(model.ref)$adj.r.squared
    
    res[1,"RSE.refchems"] <- rse.refchems
    res[1,"R2.refchems"] <- r2.refchems
    do.print <- F
    if(do.print) print(summary(model))
    if(sum(amask)>1) {
      xp <- as.matrix(mat.med[,amask==1]) %*% coefs
      xp.min <- as.matrix(mat.min[,amask==1]) %*% coefs
      xp.max <- as.matrix(mat.max[,amask==1]) %*% coefs
    }
    else {
      xp <- mat.med[,amask==1] * coefs
      xp.min <- mat.min[,amask==1] * coefs
      xp.max <- mat.max[,amask==1] * coefs
    }
    xp[xp<0] <- 0
    xp.min[xp.min<0] <- 0
    xp.max[xp.max<0] <- 0
    
    yp <- auc.model[,"AUC.agonist"]
    xp2 <- xp
    yp2 <- yp
    xp2[xp2<cutoff] <- 0
    xp2[xp2>0] <- 1
    yp2[yp2<cutoff] <- 0
    yp2[yp2>0] <- 1
    a <- sum(xp2*yp2)
    b <- sum(xp2*(1-yp2))
    c <- sum((1-xp2)*yp2)
    d <- sum((1-xp2)*(1-yp2))
    txt <- TxT(a,b,c,d)
    res[1,"sens.allchem"] <- txt$sens
    res[1,"spec.allchem"] <- txt$spec
    res[1,"ba.allchem"] <- txt$ba
    ba.all <- txt$ba
    do.print <- F
    if(ba.all>bamax || sum(amask)==length(amask)) {
      do.print <- T
      bamax <- ba.all
    }

    if(do.print) {
      main <- ""
      for(k in 1:NASSAY) main <- paste(main,amask[k],sep="")
      cat(sum(amask),":",main,"\n")
      plot(yp~xp,xlim=c(0,1),ylim=c(0,1),xlab="Subset Model AUC Score",ylab="Full Model AUC Score",main=main)
      for(j in 1:nchem) {
        x <- xp[j]
        y <- yp[j]
        minval <- auc.model[j,"minbar"]
        maxval <- auc.model[j,"maxbar"]
        lines(c(x,x),c(minval,maxval))
        if(x>=minval && x<=maxval) points(x,y,bg="red",cex=1,pch=21)
        
        minval <- xp.min[j]
        maxval <- xp.max[j]
        lines(c(minval,maxval),c(y,y))
        if(y>=minval && y<=maxval) points(x,y,bg="red",cex=1,pch=21)
        
      }
      lines(c(0,1),c(0,1))
      lines(c(0,1),c(cutoff,cutoff))
      lines(c(cutoff,cutoff),c(0,1))
      text(cutoff,0.99,paste("Assays:",sum(amask)),pos=4,cex=0.9)
      text(cutoff,0.92,paste("RMSE all chems:",format(rse.allchems,digits=2)),pos=4,cex=0.9)
      text(cutoff,0.85,paste("R2 all chems:",format(r2.allchems,digits=2)),pos=4,cex=0.9)
      text(cutoff,0.78,paste("RMSE ref chems:",format(rse.refchems,digits=2)),pos=4,cex=0.9)
      text(cutoff,0.71,paste("R2 ref chems:",format(r2.refchems,digits=2)),pos=4,cex=0.9)
      text(0.5,cutoff+0.21,"Sens,Spec,BA",pos=4,cex=1)
      text(0.5,cutoff+0.14,paste("All: [",format(txt$sens,digits=2),":",format(txt$spec,digits=2),":",format(txt$ba,digits=2),"]",sep=""),pos=4,cex=1)
    }
    xp2 <- xp2[1:nchem.ref]
    yp2 <- yp2[1:nchem.ref]
    a <- sum(xp2*yp2)
    b <- sum(xp2*(1-yp2))
    c <- sum((1-xp2)*yp2)
    d <- sum((1-xp2)*(1-yp2))
    txt <- TxT(a,b,c,d)      
    res[1,"sens.refchem"] <- txt$sens
    res[1,"spec.refchem"] <- txt$spec
    res[1,"ba.refchem"] <- txt$ba
    if(do.print) text(0.5,cutoff+0.07,paste("Ref: [",format(txt$sens,digits=2),":",format(txt$spec,digits=2),":",format(txt$ba,digits=2),"]",sep=""),pos=4,cex=1)

    code.list.invitro <- REFCHEM.DATA[!is.na(REFCHEM.DATA[,"InVitro"]),"CODE"]
    yiv <- as.numeric(REFCHEM.DATA[code.list.invitro,"InVitro"])
    if(sum(amask)>1) xiv <- xp[code.list.invitro,1]
    else xiv <- xp[is.element(code.list,code.list.invitro)]
    xiv[xiv<cutoff] <- 0
    xiv[xiv>0] <- 1
    a <- sum(xiv*yiv)
    b <- sum(xiv*(1-yiv))
    c <- sum((1-xiv)*yiv)
    d <- sum((1-xiv)*(1-yiv))
    txt <- TxT(a,b,c,d)      
    res[1,"sens.refchem.invitro"] <- txt$sens
    res[1,"spec.refchem.invitro"] <- txt$spec
    res[1,"ba.refchem.invitro"] <- txt$ba
    
    code.list.invivo <- REFCHEM.DATA[!is.na(REFCHEM.DATA[,"InVivo"]),"CODE"]
    yiv <- as.numeric(REFCHEM.DATA[code.list.invivo,"InVivo"])
    if(sum(amask)>1) xiv <- xp[code.list.invivo,1]
    else xiv <- xp[is.element(code.list,code.list.invivo)]
    xiv[xiv<cutoff] <- 0
    xiv[xiv>0] <- 1
    a <- sum(xiv*yiv)
    b <- sum(xiv*(1-yiv))
    c <- sum((1-xiv)*yiv)
    d <- sum((1-xiv)*(1-yiv))
    txt <- TxT(a,b,c,d)      
    res[1,"sens.refchem.invivo"] <- txt$sens
    res[1,"spec.refchem.invivo"] <- txt$spec
    res[1,"ba.refchem.invivo"] <- txt$ba
    
    res[1,"ba.min"] <- min(res[1,"ba.allchem"],res[1,"ba.refchem"],res[1,"ba.refchem.invitro"],res[1,"ba.refchem.invivo"])
  
    res.all <- rbind(res.all,res)
    res.amask.all <- rbind(res.amask.all,res.amask)
    res.coefs.all <- rbind(res.coefs.all,res.coefs)
    if(!to.file && do.print) browser()
    
    if(i%%1000==0) {
      dtime <- as.numeric(Sys.time()-time0,units="secs")
      fcalls <- 1000*nchem
      cat("====================================================\n")
      cat("finished: ",i," [",sum(amask),"] R2 MAX",format(bamax,digits=4),"\n")
      cat("====================================================\n")
      time0 <- Sys.time()
      flush.console()
    }
  }
  SAMPLE.STATS.ALLCHEMS.1 <<- cbind(res.all,res.amask.all,res.coefs.all)
  export.sample.stats.1(chemset,cutoff)
  if(to.file) dev.off()
}
#--------------------------------------------------------------------------------------
#
# export the sample stats
#
#--------------------------------------------------------------------------------------
summarize.sample.stats <- function() {
  res <- SAMPLE.STATS.ALLCHEMS.1
  nchem <- dim(res)[1]
  for(i in 1:nchem) {
    res[i,"ba.min"] <- min(res[i,"ba.allchem"],res[i,"ba.refchem"],res[i,"ba.refchem.invitro"],res[i,"ba.refchem.invivo"])
  }
  SAMPLE.STATS.ALLCHEMS.1 <<- res
  export.sample.stats.1("allchems",0.1)
}  
#--------------------------------------------------------------------------------------
#
# export the sample stats
#
#--------------------------------------------------------------------------------------
export.sample.stats.1 <- function(chemset,cutoff) {
  print.current.function ()
  fname <- paste("../output/sample_stats_1_",chemset,"_",cutoff,".xlsx",sep="")
  if(chemset=="allchems") write.xlsx(SAMPLE.STATS.ALLCHEMS.1,file=fname, row.names=F)
  if(chemset=="refchems") write.xlsx(SAMPLE.STATS.REFCHEMS.1,file=fname, row.names=F)
}
#--------------------------------------------------------------------------------------
#
# plot the stats - produces Figure 2 of the manuscript
#
#--------------------------------------------------------------------------------------
plot.stats.1 <- function(to.file=F) {
  print.current.function ()
  if(to.file) {
    fname <- paste("../plots/Fig 2 model_model_stats.pdf",sep="")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,1))
  data <- SAMPLE.STATS.ALLCHEMS.1
  
  groups <- data[,"Assays"]
  ba <- data[,"ba.allchem"]
  boxplot(ba~groups,xlab="Assays Included",ylab="Balanced Accuracy",cex.axis=0.9,cex.lab=1.2,ylim=c(0.7,1),main="All Chemicals")
  lines(c(0,100),c(0.9638,0.9638))
  text(15,0.75,"A",cex=2,pos=4)
  #lines(c(0,100),c(1,1))
  
  groups <- data[,"Assays"]
  ba <- data[,"ba.refchem"]
  boxplot(ba~groups,xlab="Assays Included",ylab="Balanced Accuracy",cex.axis=0.9,cex.lab=1.2,ylim=c(0.7,1),main="Reference Chemicals")
  lines(c(0,100),c(0.9638,0.9638))
  text(15,0.75,"B",cex=2,pos=4)
  #lines(c(0,100),c(1,1))
  if(!to.file) browser()
  
  groups <- data[,"Assays"]
  ba <- data[,"RSE.allchems"]
  boxplot(ba~groups,xlab="Assays Included",ylab="RSE",cex.axis=0.9,cex.lab=1.2,ylim=c(0,0.25),main="All Chemicals")
  lines(c(0,100),c(0.95,0.95))
  
  groups <- data[,"Assays"]
  ba <- data[,"RSE.refchems"]
  boxplot(ba~groups,xlab="Assays Included",ylab="RSE",cex.axis=0.9,cex.lab=1.2,ylim=c(0,0.25),main="Reference Chemicals")
  lines(c(0,100),c(0.95,0.95))
  
  if(to.file) dev.off()
  else browser()
}
#--------------------------------------------------------------------------------------
#
# add an index column to SAMPLE.STATS.ALLCHEMS.1
#
#--------------------------------------------------------------------------------------
index.sample.stats.1 <- function(do.prep=F,to.file=F) {
  print.current.function ()
  
  if(do.prep) {
    temp <- SAMPLE.STATS.ALLCHEMS.1[,ASSAY.LIST]
    index <- temp[,2]
    for(i in 2:NASSAY) index <- paste(index,temp[,i],sep="")
    all.data <- cbind(SAMPLE.STATS.ALLCHEMS.1,index)
    names(all.data)[dim(all.data)[2]] <- "index"
    ALL.DATA <<- all.data
  }
  best <- NULL
  for(i in 1:16) {
    temp <- ALL.DATA[ALL.DATA[,"Assays"]==i,]
    index <- which.max(temp[,"ba.min"])
    best <- rbind(best,temp[index,])
  }
  
  BEST.MODEL <<- best
  file <- "../output/best_model.xlsx"
  write.xlsx(BEST.MODEL,file)
  
  if(to.file) {
    fname <- paste("../plots/Fig 3 best_model.pdf",sep="")
    pdf(file=fname,width=8,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(5,4,4,4))

  x <- seq(from=1,to=16)
  x1 <- x+0.5
  y1 <- BEST.MODEL[,"ba.allchem"]
  y1a <- BEST.MODEL[,"sens.allchem"]
  y1b <- BEST.MODEL[,"spec.allchem"]
  y2 <- BEST.MODEL[,"ba.refchem"]
  y3 <- BEST.MODEL[,"ba.refchem.invitro"]
  y4 <- BEST.MODEL[,"ba.refchem.invivo"]
  plot(0~0,type="n",xlim=c(1,16),ylim=c(0.75,1),xlab="Assays",ylab="Balanced Accuracy",cex.lab=1.5,cex.axis=0.01,xaxs=NULL)
  axis(side=1,at=seq(from=1,to=16),cex.axis=1.2)
  axis(side=2,at=seq(from=0.75,to=1.0,by=0.05),cex.axis=1.5)
  lines(c(0,100),c(1.00,1.00),col="gray")
  lines(c(0,100),c(0.95,0.95),col="gray")
  lines(c(0,100),c(0.90,0.90),col="gray")
  lines(c(0,100),c(0.85,0.85),col="gray")
  points(y1~x,cex=1.5,pch=21,bg="black")
  for(i in 1:16) lines(c(i,i),c(y1a[i],y1b[i]))
  points(y1a~x,cex=1.5,pch=25,bg="black")
  points(y1b~x,cex=1.5,pch=24,bg="black")
  points(y2~x1,cex=1.5,pch=21,bg="red")
  points(y3~x1,cex=1.5,pch=22,bg="red")
  points(y4~x1,cex=1.5,pch=23,bg="red")
  

  points(4,0.825,cex=1.5,pch=21,bg="black")
  points(4,0.81,cex=1.5,pch=25,bg="black")
  points(4,0.795,cex=1.5,pch=24,bg="black")
  points(4,0.78,cex=1.5,pch=21,bg="red")
  points(4,0.765,cex=1.5,pch=22,bg="red")
  points(4,0.75,cex=1.5,pch=23,bg="red")

  text(4.5,0.825,"BA (best all chemical model vs. full model)",cex=1.1,pos=4)
  text(4.5,0.81,"Sensitivity (best all chemical model vs. full model)",cex=1.1,pos=4)
  text(4.5,0.795,"Specificity (best all chemical model vs. full model)",cex=1.1,pos=4)
  text(4.5,0.78,"BA (in vitro reference chemicals vs. full model)",cex=1.1,pos=4)
  text(4.5,0.765,"BA (in vitro reference chemicals vs. literature)",cex=1.1,pos=4)
  text(4.5,0.75,"BA (in vivo reference chemicals vs. literature)",cex=1.1,pos=4)
  
  if(to.file) dev.off()
  else browser()
}
#--------------------------------------------------------------------------------------
#
# export the ALL.DATA structure with performance above some threshold
#
#--------------------------------------------------------------------------------------
export.all.data <- function(to.file=F,threshold=0.92) {
  print.current.function ()
  temp <- ALL.DATA[ALL.DATA[,"ba.min"]>threshold,]
  file <- paste("../output/S2 all_data_",threshold,".xlsx",sep="")
  write.xlsx(temp,file)
  if(to.file) {
    fname <- paste("../plots/all_data_",threshold,".pdf",sep="")
    pdf(file=fname,width=6,height=6,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(5,4,4,4))
  hist(temp[,"Assays"],xlab="Assays",ylab="Models",main=paste("Threshold:",threshold),cex.lab=1.2,cex.axis=1.2)
  atemp <- t(temp[,ASSAY.LIST])
  rs <- rowSums(atemp)
  rs <- rs/dim(temp)[1]
  rs <- as.data.frame(rs)
  rs <- cbind(ASSAY.LIST,rs)
  names(rs) <- c("Assay","Frequency")
  file <- file <- paste("../output/all_data_assay_freq_",threshold,".xlsx",sep="")
  write.xlsx(rs,file,rowNames=F)
  if(to.file) dev.off()
  else browser()
}
#--------------------------------------------------------------------------------------
#
# get the list of chemicals that are misclassified in the best models
#
#--------------------------------------------------------------------------------------
count.misclassified.chems <- function(to.file=F,cutoff=0.1) {
  print.current.function ()
  
  data <- as.matrix(MAT.MED)
  
  for(i in 1:dim(data)[2]) {
    x <- data[,i]
    x <- x[x>0]
    med <- median(x)
    data[,i] <- data[,i]/med
  }
  data <- data/max(data)

  code.list <- FULL.MODEL[,"CODE"]
  auc <- FULL.MODEL[,"AUC.agonist"]
  rs <- rowMax(data)
  if(to.file) {
    fname <- paste("../plots/misclassified_chems_",cutoff,".pdf",sep="")
    pdf(file=fname,width=6,height=6,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(5,4,4,4))
  plot(auc~rs,xlab="Maximum Assay AUC",ylab="Full Model AUC",cex.lab=1.2,cex.axis=1.2)
  lines(c(0,1),c(0,1))
  z <- rs-auc
  res <- cbind(FULL.MODEL,z)
  names(res)[dim(res)[2]] <- "MaxAssayAUC.minus.FullModelAUC"
  file <- "../output/misclassified_chems.xlsx"
  write.xlsx(res,file)
  if(to.file) dev.off()
  else browser()
}
#--------------------------------------------------------------------------------------
#
# generate a heat map of the misclassified chemicals
#
#--------------------------------------------------------------------------------------
misclass.hm <- function(to.file=F,cutoff=0.35) {
  print.current.function ()
  if(to.file) {
    fname <- paste("../plots/Fig 4 misclass_hm.pdf",sep="")
    pdf(file=fname,width=8,height=8,pointsize=10,bg="white",paper="letter",pagecentre=T)
  }
  file <- "../output/misclassified_chems.xlsx"
  res <- read.xlsx(file)
  res <- res[res["MaxAssayAUC.minus.FullModelAUC"]>cutoff,]
  code.list <- res[,"CODE"]
  name.list <- res[,"Name"]
  
  data <- as.matrix(MAT.MED)
  
  for(i in 1:dim(data)[2]) {
    x <- data[,i]
    x <- x[x>0]
    med <- median(x)
    data[,i] <- data[,i]/med
  }
  data <- data/max(data)
  data <- data[code.list,]

  heatmap(data,margins=c(15,15),scale="none",
          xlab="",ylab="",cexCol=0.8,cexRow=0.8,col=brewer.pal(9,"Reds"),
          hclustfun=function(x) hclust(d=dist(x),method="ward.D"),labRow=name.list,
          main="")
  
  if(to.file) dev.off()
  else browser()
}

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################