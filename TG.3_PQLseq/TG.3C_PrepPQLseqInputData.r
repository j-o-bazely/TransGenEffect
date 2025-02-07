####                    TG.3C_PrepPQLseqInputData                                 ####
 ## Created by: Charley                                                            ##
 ## Date: 13.12.2024    
 ##
 ## Adapted by: James
 ## Date: 28.01.24
 ##
 ## Context: PRepares the methylKit uniteCov object                                ##  
 ##          for PQLseq requirements                                               ##  
 ##          Inc. split by chromosome: Makes the process MUCH faster, much faster  ##
####                                                                              ####

#### Prep Environment ####
## Load packages ##
library(methylKit)
library(GenomicRanges)
#library(genomation)
library(tidyverse)
library(readxl)

## Set directories ##
DIR_UniteCov <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.1_data"
DIR_PQL <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.3_data"
DIR_CHROM <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.3_data"

## Load functions
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

#### TP1 ####
#### DS  ####
TP1DS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'DS', TP = 'TP1')
#### MD  ####
TP1MD <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MD', TP = 'TP1')
#### MS  ####
TP1MS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MS', TP = 'TP1')

#### TP3 ####
#### DS  ####
TP3DS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'DS', TP = 'TP3')
#### MD  ####
TP3MD <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MD', TP = 'TP3')
#### MS  ####
TP3MS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MS', TP = 'TP3')


#### TP6 ####
#### DS  ####
TP6DS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'DS', TP = 'TP6')
#### MD  ####
TP6MD <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MD', TP = 'TP6')
#### MS  ####
TP6MS <- PrepInputObjs(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT = 'MS', TP = 'TP6')
