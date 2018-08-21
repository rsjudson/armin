#--------------------------------------------------------------------------------------
#'
#' prep the T matrix with a variable number of assays
#'
#' Pathway-specific - update pathway-specific TMAT
#'
#' TMAT is a matrix whose elements match the connectivity in the pathway drawing
#--------------------------------------------------------------------------------------
TmatVa <- function() {
  #printCurrentFunction()
  temp <- matrix(nrow=NASSAY,ncol=NRECEPTOR)
  temp[] <- 0
  ##CHANGED
  #Agonist R1
  temp[,1] <-  c(1,1,1,1,1,1,1,1,1,0,0,1,0)
  #Antagonist R2
  temp[,2] <-  c(1,1,1,1,1,0,0,0,0,1,1,0,1)
  #R3 through R7
  temp[,3] <-  c(1,1,1,0,0,0,0,0,0,0,0,0,0)
  temp[,4] <-  c(0,0,0,1,1,0,0,0,0,0,0,0,0)
  temp[,5] <-  c(0,0,0,0,0,1,0,0,0,0,0,0,0)
  temp[,6] <-  c(0,0,0,0,0,0,1,1,1,0,0,0,0)
  temp[,7] <-  c(0,0,0,0,0,0,0,0,0,1,1,0,0)
  #Assay by assay interference
  # NVS 
  temp[,8] <- c(1,0,0,0,0,0,0,0,0,0,0,0,0)
  temp[,9] <- c(0,1,0,0,0,0,0,0,0,0,0,0,0)
  temp[,10] <- c(0,0,1,0,0,0,0,0,0,0,0,0,0)
  # OT PC
  temp[,11] <- c(0,0,0,1,0,0,0,0,0,0,0,0,0)
  temp[,12] <- c(0,0,0,0,1,0,0,0,0,0,0,0,0)
  # ATG not needed (same as R5)
  # OT ARE
  temp[,13] <- c(0,0,0,0,0,0,1,0,0,0,0,0,0)
  # NCGC Agonist
  temp[,14] <- c(0,0,0,0,0,0,0,1,0,0,0,0,0)
  temp[,15] <- c(0,0,0,0,0,0,0,0,1,0,0,0,0)
  # NCGC Antagonist
  temp[,16] <- c(0,0,0,0,0,0,0,0,0,1,0,0,0)
  temp[,17] <- c(0,0,0,0,0,0,0,0,0,0,1,0,0)

  # ACEA Agonist
  temp[,18] <- c(0,0,0,0,0,0,0,0,0,1,0,1,0)
  
  # ACEA Antagonist
  temp[,19] <- c(0,0,0,0,0,0,0,0,0,1,0,0,1)

  temp <- temp[AMASK==1,]
  rmask <- colSums(temp)
  rmask[rmask>0] <- 1
  temp <- temp[,rmask==1]
  TMAT <<- temp
  RMASK <<- rmask
}
