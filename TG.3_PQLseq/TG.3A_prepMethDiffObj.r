### TG.3A_prepMethDiffObj.r ###
## Created: February 7th 2025
## By : James Bazely
## Context: This script generates diffmeth objects between each intersection,
##          within each timepoint


## load function/s
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

## Set directories ##
DIR_DATA <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage"
DIR_PQL <- file.path(DIR_DATA, "TG_analysisData/TG.3_data")


prepMethDiffObj <- function(TP, INTER){
  ## Load object
  intersectObj <- readRDS(file.path(paste(DIR_DATA, '/TG_analysisData/TG.1_data/', TP, '/TG_', TP, '_', INTER, 'methObj75pc.RDS', sep = '')))
  
  ## Generate methDiff object
  diffMethObj <- calculateDiffMeth(intersectObj, mc.cores = 10)
  
  ## Save methDiff object
  saveRDS(diffMethObj, file.path(paste(DIR_DATA, '/TG_analysisData/TG.3_data/', TP, '/TG_', TP, '_', INTER, 'diffMethObj75pc.RDS', sep = '')))
}

timepoints <- c('TP1', 'TP3', 'TP6')
intersections <- c('DS', 'MD', 'MS')

for(i in timepoints){
  for(j in intersections){
    print(paste('Generating diffMeth object for: ', i, '_', j, sep = ''))
    prepMethDiffObj(TP = i, INTER = j)
    gc()
  }
}










