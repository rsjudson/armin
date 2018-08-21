#--------------------------------------------------------------------------------------
#'
#' calculate the calcibration curve between AC50 and AUC
#'
#--------------------------------------------------------------------------------------
calibrationCurve <- function(to.file=F) {
  printCurrentFunction()
  ac50s <- c(0.00001,0.00002,0.00005,
             0.0001,0.0002,0.0005,
             0.001,0.002,0.005,
             0.01,0.02,0.05,
             0.1,0.2,0.5,
             1,2,5,
             10,20,50,
             100,200,500,
             1000,2000,5000,10000,100000)
  nac50 <- length(ac50s)
  aucs <- ac50s
  aucs[] <- 0
  for(i in 1:nac50) {
    ac50 <- ac50s[i]
    top <- 1  ##Specific to antagonist pathway: use 0.55
    w <- 1
    cr.mat <- vector(length=length(CONCLIST),mode="numeric")
    cr.mat[] <- 0
    for(j in 1:length(CONCLIST)) {
      conc <- CONCLIST[j]
      cr.mat[j] <- top*(conc**w/(conc**w+ac50**w))
    }
    aucs[i] <- receptorScore(cr.mat,method=2,do.print=F)
    cat("ac50,auc,",ac50s[i],":",format(aucs[i],digits=2),"\n")
  }
  if(to.file) {
    fname <- "../plots/calibration_curve.pdf"
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(6,6,3,3))
  
  plot(aucs~ac50s,xlab="AC50 (uM)",ylab="AUC",xlim=c(1e-4,1e4),log="xy",cex.lab=1.2, cex.axis=1.2,type="l",ylim=c(1e-3,1.5),xaxp=c(0.0001,1000,1),yaxp=c(0.0001,1,1))
  lines(c(1e-7,1e7),c(0.001,0.001),col="gray")
  lines(c(1e-7,1e7),c(0.01,0.01),col="gray")
  lines(c(1e-7,1e7),c(0.1,0.1),col="gray")
  lines(c(1e-7,1e7),c(1,1),col="gray")
  
  lines(c(1e-4,1e-4),c(1e-6,10),col="gray")
  lines(c(1e-3,1e-3),c(1e-6,10),col="gray")
  lines(c(1e-2,1e-2),c(1e-6,10),col="gray")
  lines(c(1e-1,1e-1),c(1e-6,10),col="gray")
  lines(c(1,1),c(1e-6,10),col="gray")
  lines(c(10,10),c(1e-6,10),col="gray")
  lines(c(100,100),c(1e-6,10),col="gray")
  lines(c(1000,1000),c(1e-6,10),col="gray")
  lines(c(10000,10000),c(1e-6,10),col="gray")
  
  if(to.file) dev.off()
  else browser()
  
}
