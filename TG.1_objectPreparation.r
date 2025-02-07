####                    TG.1_ObjectPrep.R                                   ####
## Created: December 2nd 2024
## By : James Bazely
## Context: This script prepares 100% coverage objects for each key analysis group; 
##          mothers, hatchlings, deep hatchlings and shallow hatchlings.
##          It also prepares these objects split into the three time points used 
##          within the split-clutch experiment, representing three different temperatures
##          across the nesting season
##          Next, it creates the necessary GRanges objects to interrogate each intersection
##          at each of the three time points. 

### Note: - SLM samples can't be included in time point analysis (samples are from 2022)



#### Timepoint 1 objects                                                    ####
### Env Prep ###
## Load libraries ##
library(tidyverse)
library(methylKit)

## Set number of cores ##
nslots = 1

## Set paths ##
dataDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage'

## load function/s
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

## Load metadata
metadata <- read.csv(file.path(dataDir, 'TG_metadata/Metadata_TSD_All_Hatchlings.csv')) # Featuring all-new metadata!

      #### MethylKit Object Prep ####
### Mother methylKit ###
TP1_methObj_mother <- prepMethObj75pc(metadata,
                                         group = 'mother',
                                         tp = '1')
## Save ##
saveRDS(TP1_methObj_mother, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_motherMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Hatchling methylKit ###
TP1_methObj_hatch <- prepMethObj75pc(metadata,
                                        group = 'hatchling',
                                        tp = '1')
## Save ##
saveRDS(TP1_methObj_hatch, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_hatchMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Deep Hatchling intersection methylKit ###
TP1_methObj_deep <- prepMethObj75pc(metadata, 
                                       group = 'deep',
                                       tp = '1')
## Save ##
saveRDS(TP1_methObj_deep, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_deepMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Shallow hatchling intersection methylKit ###
TP1_methObj_shallow <- prepMethObj75pc(metadata, 
                                          group = 'shallow',
                                          tp = '1')
## Save ##
saveRDS(TP1_methObj_shallow, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_shallowMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Deep / Shallow intersection methylKit
TP1_methObj_DS <- PrepIntersectObj75pc(metadata,
                                   group = 'DS',
                                   tp = '1')
## Save ##
saveRDS(TP1_methObj_DS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_DSmethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Deep / Shallow intersection methylKit ###
TP1_MethObj_MS <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MS',
                                        tp = '1')
## Save ##
saveRDS(TP1_MethObj_MS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_MSmethObj75 pc.RDS'))

## Clear free RAM ##
gc()


### Mother / Deep intersection methylKit ###
TP1_methObj_MD <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MD',
                                        tp = '1')
## Save ##
saveRDS(TP1_methObj_MD, file.path(dataDir, 'TG_analysisData/TG.1_data/TP1/TG_TP1_MDmethObj75pc.RDS'))

## Clear free RAM ##
gc()


            #### GRanges Object Prep ####
prepGlobMethObj75pc(tp = 'TP1')


#### Timepoint 3 objects                                                    ####
### Env Prep ###
## Load libraries ##
library(tidyverse)
library(methylKit)

## Set number of cores ##
nslots = 4

## Set paths ##
dataDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage'

## load function/s
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

## Load metadata
metadata <- read.csv(file.path(dataDir, 'TG_metadata/Metadata_TSD_All_Hatchlings.csv')) # Featuring all-new metadata!

      #### MethylKit Object Prep ####
### Mother methylKit ###
TP3_methObj_mother <- prepMethObj75pc(metadata,
                                       group = 'mother',
                                       tp = '3')

## Save object ##
saveRDS(TP3_methObj_mother, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_motherMethObj75pc.RDS'))

## Clear free RAM ##
gc()

### Hatchling methylKit ###
TP3_methObj_hatch <- prepMethObj75pc(metadata,
                                      group = 'hatchling',
                                      tp = '3')

## Save object ##
saveRDS(TP3_methObj_hatch, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_hatchMethObj75pc.RDS'))

## Clear free RAM ##
gc()

### Deep Hatchling methylKit ###
TP3_methObj_deep <- prepMethObj75pc(metadata, 
                                     group = 'deep',
                                     tp = '3')

## Save object ##
saveRDS(TP3_methObj_deep, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_deepMethObj75pc.RDS'))

# Clear free RAM ##
gc()

### Shallow hatchling methylKit ###
TP3_MethObj_shallow <- prepMethObj75pc(metadata, 
                                          group = 'shallow',
                                          tp = '3')

## Save object ##
saveRDS(TP3_MethObj_shallow, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_shallowMethObj75pc.RDS'))

## Clear free RAM
gc()


### Deep / Shallow intersection methylKit ###
TP3_methObj_DS <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'DS',
                                        tp = '3')

## Save object ##
saveRDS(TP3_methObj_DS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_DSmethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Mother / Deep intersection methylKit ###
TP3_methObj_MD <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MD',
                                        tp = '3')

## Save object ##
saveRDS(TP3_methObj_MD, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_MDmethObj75pc.RDS'))

## Clear free RAM ##
gc()

### Mother / Shallow intersection methylKit ###
TP3_methObj_MS <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MS',
                                        tp = '3')


### Saving objects
saveRDS(TP3_methObj_MS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP3/TG_TP3_MSmethObj75pc.RDS'))

## Clear free RAM ##
gc()

            #### GRanges Object Prep ####
prepGlobMethObj75pc(tp = 'TP3')

#### Timepoint 6 objects                                                    ####
### Env Prep ###
## Load libraries ##
library(tidyverse)
library(methylKit)

## Set number of cores ##
nslots = 4

## Set paths ##
dataDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage'

## load function/s
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

## Load metadata
metadata <- read.csv(file.path(dataDir, 'TG_metadata/Metadata_TSD_All_Hatchlings.csv')) # Featuring all-new metadata!


      #### MethylKit Object Prep ####
### Mother methylKit ###
TP6_methObj_mother <- prepMethObj75pc(metadata,
                                       group = 'mother',
                                       tp = '6')

saveRDS(TP6_methObj_mother, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_motherMethObj75pc.RDS'))
 
## Clear free RAM ##
gc()


### Hatchling methylKit ###
TP6_methObj_hatch <- prepMethObj75pc(metadata, 
                                     group = 'hatchling',
                                     tp = '6'
)

saveRDS(TP6_methObj_hatch, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_hatchMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Deep Hatchling methylKit ###
TP6_methObj_deep <- prepMethObj75pc(metadata, 
                                     group = 'deep',
                                     tp = '6'
)

saveRDS(TP6_methObj_deep, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_deepMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Shallow hatchling methylKit ###
TP6_methObj_shallow <- prepMethObj75pc(metadata, 
                                          group = 'shallow',
                                          tp = '6')

saveRDS(TP6_methObj_shallow, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_shallowMethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Deep / Shallow intersection methylKit ###
TP6_methObj_DS <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'DS',
                                        tp = '6')

saveRDS(TP6_methObj_DS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_DSmethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Mother / Deep intersection methylKit ###
TP6_methObj_MD <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MD',
                                        tp = '6')


saveRDS(TP6_methObj_MD, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_MDmethObj75pc.RDS'))

## Clear free RAM ##
gc()


### Mother / Shallow intersection methylKit ###
TP6_methObj_MS <- PrepIntersectObj75pc(metadata = metadata,
                                        group = 'MS',
                                        tp = '6')

saveRDS(TP6_methObj_MS, file.path(dataDir, 'TG_analysisData/TG.1_data/TP6/TG_TP6_MSmethObj75pc.RDS'))

## Clear free RAM ##
gc()

            #### GRanges Object Prep ####
prepGlobMethObj75pc(tp = 'TP6')

#### Combined timepoint objects                                             ####
#### MethylKit Object Prep ####
### Mother methylKit ###
combo_methObj_mother <- prepMethObj75pc(metadata,
                                         group = 'mother',
                                         tp = 'all_three!!!')

### Hatchling methylKit ###
combo_methObj_hatch <- prepMethObj75pc(metadata,
                                        group = 'hatchling', 
                                        tp = 'all_three!!!')

### Deep Hatchling methylKit ###
combo_methObj_deep <- prepMethObj75pc(metadata, 
                                       group = 'deep',
                                       tp = 'all_three!!!'
)

### Shallow hatchling methylKit ###
combo_methObj_shallow <- prepMethObj75pc(metadata, 
                                          group = 'shallow',
                                          tp = 'all_three!!!')

### Saving objects
saveRDS(combo_methObj_mother, file.path(dataDir, 'TG_analysisData/TG.1_data/TG_Combo/TG_combo_motherMethObj75pc.RDS'))
saveRDS(combo_methObj_hatch, file.path(dataDir, 'TG_analysisData/TG.1_data/TG_Combo/TG_combo_hatchMethObj75pc.RDS'))
saveRDS(combo_methObj_deep, file.path(dataDir, 'TG_analysisData/TG.1_data/TG_Combo/TG_combo_deepMethObj75pc.RDS'))
saveRDS(combo_methObj_shallow, file.path(dataDir, 'TG_analysisData/TG.1_data/TG_Combo/TG_combo_shallowMethObj75pc.RDS'))






####### GETTING NUMBERS #######
TP3path <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.1_data/TG_TP3'
TG_TP3_M <- makeGrangeFromMethObj(readRDS(file.path(TP3path, 'TG_TP3_motherMethObj75pc.RDS')))
TG_TP3_H <- makeGrangeFromMethObj(readRDS(file.path(TP3path, 'TG_TP3_hatchMethObj75pc.RDS')))
TG_TP3_S <- makeGrangeFromMethObj(readRDS(file.path(TP3path, 'TG_TP3_shallowMethObj75pc.RDS')))
TG_TP3_D <- makeGrangeFromMethObj(readRDS(file.path(TP3path, 'TG_TP3_deepMethObj75pc.RDS')))
TG_TP3_MH <- readRDS(file.path(TP3path, 'TG_TP3_MHgrange.RDS'))
TG_TP3_MS <- readRDS(file.path(TP3path, 'TG_TP3_MSgrange.RDS'))
TG_TP3_MD <- readRDS(file.path(TP3path, 'TG_TP3_MDgrange.RDS'))
TG_TP3_SD <- readRDS(file.path(TP3path, 'TG_TP3_SDgrange.RDS'))
TG_TP3_SDM <- readRDS(file.path(TP3path, 'TG_TP3_SDMgrange.RDS'))
TG_TP3_MDexc <- readRDS(file.path(TP3path, 'TG_TP3_MDgrange_exclusive.RDS'))
TG_TP3_SDexc <- readRDS(file.path(TP3path, 'TG_TP3_SDgrange_exclusive.RDS'))
TG_TP3_MSexc <- readRDS(file.path(TP3path, 'TG_TP3_MSgrange_exclusive.RDS'))

###

TP6path <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.1_data/TG_TP6'
TG_TP6_M <- makeGrangeFromMethObj(readRDS(file.path(TP6path, 'TG_TP6_motherMethObj75pc.RDS')))
TG_TP6_H <- makeGrangeFromMethObj(readRDS(file.path(TP6path, 'TG_TP6_hatchMethObj75pc.RDS')))
TG_TP6_S <- makeGrangeFromMethObj(readRDS(file.path(TP6path, 'TG_TP6_shallowMethObj75pc.RDS')))
TG_TP6_D <- makeGrangeFromMethObj(readRDS(file.path(TP6path, 'TG_TP6_deepMethObj75pc.RDS')))
TG_TP6_MH <- readRDS(file.path(TP6path, 'TG_TP6_MHgrange.RDS'))
TG_TP6_MS <- readRDS(file.path(TP6path, 'TG_TP6_MSgrange.RDS'))
TG_TP6_MD <- readRDS(file.path(TP6path, 'TG_TP6_MDgrange.RDS'))
TG_TP6_SD <- readRDS(file.path(TP6path, 'TG_TP6_SDgrange.RDS'))
TG_TP6_SDM <- readRDS(file.path(TP6path, 'TG_TP6_SDMgrange.RDS'))
TG_TP6_MDexc <- readRDS(file.path(TP6path, 'TG_TP6_MDgrange_exclusive.RDS'))
TG_TP6_SDexc <- readRDS(file.path(TP6path, 'TG_TP6_SDgrange_exclusive.RDS'))
TG_TP6_MSexc <- readRDS(file.path(TP6path, 'TG_TP6_MSgrange_exclusive.RDS'))
