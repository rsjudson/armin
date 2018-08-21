#--------------------------------------------------------------------------------------
#'
#' Calculate the least squares solution for one chemical
#'
#' Pathway-specific - update pathway-specific assay behavior
#'
#'-------------------------------------------------------------------------------------
LsLegend <- function(to.file=F) {
  printCurrentFunction()
  if(to.file) {
    fname <- "../../manuscript/model_legend.pdf"
    pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  
  ytop <- 2.9
  plot(1~1,type="n",col.axis="white",tcl=0.01,cex.axis=0.1,cex.lab=0.1,xlim=c(0,5),ylim=c(0.8,ytop),xlab="",ylab="",lwd=3,main="")
  text(x=1,y=ytop,labels="Assay Legend",pos=4,cex=1.2)
  dy <- 0.085
  ytop <- ytop - 0.05
  ##CHANGED  
  anames <- vector(length=NASSAY,mode="character")
  anames[1] <- "A1: human AR cell-free radioligand binding (NVS)"
  anames[2] <- "A2: chimp AR cell-free radioligand binding (NVS)"
  anames[3] <- "A3: rat AR cell-free radioligand binding (NVS)"
  anames[4] <- "A4: AR-SRC protein complementation/FRET 8 hr (OT)"
  anames[5] <- "A5: AR-SRC protein complementation/FRET 16 hr (OT)"
  anames[6] <- "A6: AR-TRANS reporter gene (ATG)"
  anames[7] <- "A7: AR-ARE luciferase agonist reporter gene 24 hr (OT)"
  anames[8] <- "A8: AR beta-lactamase agonist reporter gene (Tox21)"
  anames[9] <- "A9: AR luciferase-MDAKB2 agonist reporter gene (Tox21)"
  anames[10] <- "A10: AR beta-lactamase antagonist reporter gene (Tox21)"
  anames[11] <- "A11: AR luciferase-MDAKB2 antagonist reporter gene (Tox21)"
  anames[12] <- "A12: AR Proliferation (ACEA)"
  anames[13] <- "A13: AR Proliferation suppression (ACEA)"
  ##CHANGED   	
  cols <- c("black","black","black","green","green","hotpink","cyan","cyan","cyan","orange","orange","gray","gray")
  ltys <- c(1,2,3, 1,2, 1, 1,2,3, 1,2,1,2)
  
  for(i in 1:NASSAY) {
    lines(x=c(0.01,1),y=c(ytop-dy*i,ytop-dy*i),col=cols[i],lwd=3,lty=ltys[i])
    text(x=1,y=ytop-dy*i,labels=anames[i],pos=4,cex=0.9)
  }
  ##CHANGED
  ltys <- c( 1,     1,     1,       1,      1,        1,     1,       1,     2)
  cols <- c("blue","red", "black", "green", "hotpink","cyan","orange","gray","gray")
  rnames <- vector(length=NRECEPTOR,mode="character")
  rnames[1] <- "R1: AR Agonist Model"
  rnames[2] <- "R2: AR Antagonist Model"
  rnames[3] <- "R3: Interference: cell-free radioligand binding (NVS)"
  rnames[4] <- "R4: Interference: protein complementation (PCA)/FRET (OT)"
  rnames[5] <- "R5: Interference: RNA reporter gene agonist (ATG)"
  rnames[6] <- "R6: Interference: protein reporter gene agonist (OT/Tox21)"
  rnames[7] <- "R7: Interference: protein reporter antagonist (Tox21)"
  rnames[8] <- "R18: Interference: proliferation agonist (ACEA)"
  rnames[9] <- "R19: Interference: proliferation antagonist (ACEA)"
  ytop <- 1.6
  text(x=1,y=ytop,labels="Receptor Legend",pos=4,cex=1.2)
  ytop <- ytop-0.05
  for(i in 1:NRECEPTOR0) {
    lines(x=c(0.01,1),y=c(ytop-dy*i,ytop-dy*i),col=cols[i],lwd=3,lty=ltys[i])
    text(x=1,y=ytop-dy*i,rnames[i],pos=4,cex=0.9)
  }
  if(to.file) dev.off()
}
