#--------------------------------------------------------------------------------------
#
# lane plot of AUC for reference chemicals
#
#===============================================================================
# Pathway-specific - update pathway-specific assay behavior
#===============================================================================
#--------------------------------------------------------------------------------------
refchemLaneplot <- function(to.file=F,zcut=F) {
  if(to.file) {
    fname <- paste("../plots/refchem_laneplot_",ALPHA,"_",RS0,".pdf",sep="")
    pdf(file=fname,width=7,height=6,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  file <- paste0("../output/refchem_AUC_",ALPHA,"_",RS0,".xlsx")
  rdata <- read.xlsx(file)
  nchem <- dim(rdata)[1]
  dmat <- matrix(nrow=nchem,ncol=NRECEPTOR0)
  counter <- 0
  dmat <- rdata[,4:dim(rdata)[2]]
  dmat <- dmat[,1:NRECEPTOR0]
  
  dmat[dmat<0.01] <- 0.01
  
  pot.agon <- REFCHEMS[,"agonist_potency"]
  pot.agon <- pot.agon[!is.na(pot.agon)]
  nagon <- length(pot.agon)
  
  plot(1~1,log="x",xlim=c(0.01,200),ylim=c(0,nagon),xlab="AR Pathway Model (R1)",ylab="",type="n",main="Agonist",yaxt="n",xaxp=c(0.01,1,1),cex.axis=0.8)
  x <- 1
  lines(c(x,x),c(0,nagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nagon),col="gray"); x <- x/10
  
  order.list <- sort(dmat[,"R1"],index.return=T)$ix
  counter <- 0
  for(i in 1:nchem) {
    cname <- REFCHEMS[order.list[i],"name"]
    pot.agon <- REFCHEMS[order.list[i],"agonist_potency"]
    v1 <- 1e-6
    v2 <- 1e-6
    if(!is.na(pot.agon)) {
      counter <- counter+1
      v1 <- dmat[order.list[i],"R1"]
      v2 <- max(dmat[order.list[i],])
      col="green"
      if(pot.agon=="Negative") col <- "red"
      # if(pot.agon=="Inactive") col <- "red"
      points(v1,counter,pch=21,bg=col,fg="black",cex=1.5)
      ##Comment/uncomment next 5 lines to control pseudo-receptor plotting 
      #     		if(v2>v1*2) {
      # 	    		points(v2,counter,pch=4)
      #     			ival <-  which.max( is.element(dmat[order.list[i],],max(dmat[order.list[i],])) )
      # 				rval <- names(dmat)[ival]
      #     			text(v2,counter,rval,pos=2,cex=0.8)
      #     		}
      text(1.5,counter,pot.agon,pos=4,cex=0.8)
      text(8,counter,cname,pos=4,cex=0.8)
    }
    
  }
  if(!to.file) browser()
  
  pot.antagon <- REFCHEMS[,"antagonist_potency"]
  pot.antagon <- pot.antagon[!is.na(pot.antagon)]
  nantagon <- length(pot.antagon)
  
  delta <- nantagon/nagon
  plot(1~1,log="x",xlim=c(0.01,200),ylim=c(0,nantagon),xlab="AR Pathway Model (R2)",ylab="",type="n",main="Antagonist",yaxt="n",xaxp=c(0.01,1,1),cex.axis=0.8)
  x <- 1
  lines(c(x,x),c(0,nantagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nantagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nantagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nantagon),col="gray"); x <- x/10
  lines(c(x,x),c(0,nantagon),col="gray"); x <- x/10
  order.list <- sort(dmat[,"R1"],index.return=T)$ix
  dmat <- dmat[order.list,]
  refchems <- REFCHEMS[order.list,]
  order.list <- sort(dmat[,"R2"],index.return=T)$ix
  
  lines(c(1e-6,1e6),c(nantagon+1,nantagon+1))
  counter <- 0
  for(i in 1:nchem) {
    cname <- refchems[order.list[i],"name"]
    pot.antagon <- refchems[order.list[i],"antagonist_potency"]
    v1 <- 1e-6
    v2 <- 1e-6
    if(!is.na(pot.antagon)) {
      counter <- counter+1
      v1 <- dmat[order.list[i],"R2"]
      v2 <- max(dmat[order.list[i],])
      col="green"
      if(pot.antagon=="Negative") col <- "red"
      points(v1,counter,pch=21,bg=col,fg="black",cex=1.5)
      ##Comment/uncomment next 5 lines to control pseudo-receptor plotting 
      #     		if(v2>v1*2) {
      # 	    		points(v2,counter,pch=4)
      #     			ival <-  which.max( is.element(dmat[order.list[i],],max(dmat[order.list[i],])) )
      # 				rval <- names(dmat)[ival]
      #     			text(v2,counter,rval,pos=2,cex=0.8)
      #     		}
      text(1.5,counter,pot.antagon,pos=4,cex=0.8)
      text(11,counter,cname,pos=4,cex=0.8)
    }
    
  }
  
  if(!to.file) browser()
  else dev.off()
}
