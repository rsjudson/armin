source("prepChems.R")
source("prepFiles.R")
source("prepCytotox.R")
source("prepCytotoxFilteredAC50.R")
source("prepFlags.R")

source("AC50HM.R")
source("applyFilters.R")
source("calibrationCurve.R")
source("loadData.R")
source("LsLegend.R")
source("prepZ.R")
source("TmatVa.R")

source("receptorScore.R")
source("prepCRref.R")
source("prepCRall.R")
source("prepCR.R")
source("LSref.R")
source("LSall.R")
source("LSoneVa.R")
source("interpolateConc.R")
source("AFRVa.R")
source("penalty.R")
source("transfer.R")
source("pseudoAC.R")
source("RHM.R")
source("AUCcalc.R")
source("AUCHM.R")
source("refchemLaneplot.R")
source("refchemDist.R")
source("EAScalc.R")
source("arIndex.R")
source("prepTox21Shift.R")

source("plotPenalty.R")
source("prepSupermatrix.R")
source("addSpecificityScore.R")
source("addAntagonistShift.R")
source("addAntagonistScore.R")
source("specificityScores.R")
source("allstats.R")
source("dxSpecificityZT.R")
source("dxStructureSpecificityZT.R")
source("dxPainsZtfiltered.R")
source("dxAtgCrosstalk.R")
source("dxNvsCrosstalk.R")
source("dxTox21Crosstalk.R")
source("dxFlagFilter.R")
source("dxPhyschemSpecificityZtfiltered.R")
source("receptorTree.R")
source("sigmoid.R")
source("hill.R")
source("hitClassify.R")
source("tox21AntagonistCytoFilter.R")
source("topHist.R")

source("general_utils.R")
source("db_utils.R")
source("stats_utils.R")
#===============================================================================
# Pathway-specific - start
#===============================================================================
PATHWAY <<- "AR"
NASSAY <<- 13
NRECEPTOR0 <<- 9  ##Number on pathway
NRECEPTOR <<- 19  ##NRECEPTOR + # ASSAYS (-1 bc ATG maps to one R)
AMASK <<- c(1,1,1,1,1,1,1,1,1,1,1,1,1)
#SCALE <<- 1/0.36 # value is adjusted to give AUC(agonist)=1 for positive control 
SCALE <<- 2.78  #4.173554 #previously 4.208, 5.05 - value is adjusted to give AUC(antagonist)=1 for positive control
SPECIFIC.AUC.CUTOFF.AG <<- 0.1
SPECIFIC.AUC.CUTOFF.ANT <<- 0.05
#===============================================================================
# Pathway-specific - stop
#===============================================================================

HEATMAP.CMETHOD <<- "ward.D"
NCONC <<- 14
CONCLIST <- c(0.0122,0.0244,0.0488,0.0977,0.195,0.391,0.781,1.56,3.125,6.25,12.5,25,50,100)
ALPHA <<- 0.05
WIDTH <<- 1
RS0 <<- 0.5
#--------------------------------------------------------------------------------------
#
# run all
#
#--------------------------------------------------------------------------------------
driver <- function(do.file.prep=F,
                   do.prep=F,
                   do.ref=F,
                   do.all=F,
                   do.EAS=F,
                   do.summary=F) {
  date.string <- "180821"
  if(do.file.prep) {
    prepChems(date.string)
    prepFiles(date.string,fix.nvs=F)
    tox21AntagonistCytoFilter(date.string)
    prepCytotox()
    prepCytotoxFilteredAC50()
    prepFlags(date.string)
  }
  if(do.prep) {
    loadData() 
    prepZ()
    hitClassify()
    applyFilters()
    prepTox21Shift(date.string)
    TmatVa()
    LsLegend(to.file=T)
    AC50HM(to.file=T,zcut=F)
    calibrationCurve(to.file=T)
    topHist(T)
  }
  if(do.ref) {
    prepCRref(zcut=F)
    LSref(to.file=T,zcut=F)
    RHM(mode="ref",to.file=T,zcut=F)
    AUCcalc(mode="ref",zcut=F)
    AUCHM(mode="ref",to.file=T,zcut=F,dcut=5)
    refchemLaneplot(to.file=T,zcut=F)
  }
  if(do.all) {
    prepCRall(zcut=F)
    LSall(to.file=T,zcut=F)
    AUCcalc(mode="all",zcut=F)
    AUCHM(mode="all",to.file=T,zcut=F,dcut=5)
   }
  if(do.EAS) {
    EAScalc("FP","R1")
    EAScalc("FP","R2")
  }
  if(do.summary) {
    prepSupermatrix()
    addSpecificityScore()
    addAntagonistShift()
    addAntagonistScore()
    allstats()
    #dxSpecificityZT()
    #dxStructureSpecificityZT(super=F,cutoff=0.5)
    #dxPainsZtfiltered(cutoff=0.5)
    #dxAtgCrosstalk()
    #dxNvsCrosstalk()
    #dxTox21Crosstalk()
    #dxFlagFilter()
    #dxPhyschemSpecificityZtfiltered(cutoff=0.5)
    
    #comp.lit(nmin=4)
    receptorTree()
    
  }
}
