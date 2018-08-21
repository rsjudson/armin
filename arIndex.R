#--------------------------------------------------------------------------------------
#
# give various identifiers for each receptor
#
#===============================================================================
# Pathway-specific - start
#===============================================================================
#--------------------------------------------------------------------------------------
arIndex <- function(receptor) {
  al <- ASSAY.LIST
    
  nickname <- receptor
  if(receptor=="R1") {
    nickname <- "Agonist"
    assay.list <- c(al[1],al[2],al[3],al[4],al[5],al[6],al[7],al[8],al[9],al[12])
  }
  if(receptor=="R2") {
    nickname <- "Antagonist"
    assay.list <- c(al[1],al[2],al[3],al[4],al[5],al[10],al[11],al[13])
  }
  if(receptor=="R3") {
    assay.list <- c(al[1],al[2],al[3])
  }
  if(receptor=="R4") {
    assay.list <- c(al[4],al[5])
  }
  if(receptor=="R5") {
    assay.list <- c(al[6])
  }
  if(receptor=="R6") {
    assay.list <- c(al[7],al[8],al[9])
  }
  if(receptor=="R7") {
    assay.list <- c(al[10],al[11])
  }
  if(receptor=="R8") {
    nickname <- "A1"
    assay.list <- c(al[1])
  }
  if(receptor=="R9") {
    nickname <- "A2"
    assay.list <- c(al[2])
  }
  if(receptor=="R10") {
    nickname <- "A3"
    assay.list <- c(al[3])
  }
  if(receptor=="R11") {
    nickname <- "A4"
    assay.list <- c(al[4])
  }
  if(receptor=="R12") {
    nickname <- "A5"
    assay.list <- c(al[5])
  }
  if(receptor=="R13") {
    nickname <- "A7"
    assay.list <- c(al[7])
  }
  if(receptor=="R14") {
    nickname <- "A8"
    assay.list <- c(al[8])
  }
  if(receptor=="R15") {
    nickname <- "A9"
    assay.list <- c(al[9])
  }
  if(receptor=="R16") {
    nickname <- "A10"
    assay.list <- c(al[10])
  }
  if(receptor=="R17") {
    nickname <- "A11"
    assay.list <- c(al[11])
  }
  if(receptor=="R18") {
    nickname <- "A12"
    assay.list <- c(al[12])
  }
  if(receptor=="R19") {
    nickname <- "A13"
    assay.list <- c(al[13])
  }
  return(list(nickname=nickname,assay.list=assay.list))
}
