#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 10
#$ -l h_vmem=1G
#$ -l h_rt=1:00:00
#$ -t 1-29
#$ -tc 1
##############################################################################################

# Created by Charley, 19th Mar 2024
# Run PQLSeq by chromosome in an array

##############################################################################################

### TG.3E_runPQLseq_script.sh ###

### Prep environment ###

module load R/4.2.2

SCRIPT_DIR=/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.3_PQLseq


### Call R script ###

# With 2 arguments to pass into R script:

# 1. Chromosome number (which is the job array number i.e. SGE_TASK_ID). This must be the 1st arg.
# Note we have chromosomes 0:28 but 0 can't be an array ID so the array is 1:29, 
# then in the R script the chr number is set to SGE_TASK_ID - 1

# 2. Number of cores available for script

## IMPORTANT: TIMEPOINT AND INTERSECTION VARIABLES NEED TO BE PROVIDED AFTER ${SGE_TASK_ID} ${NSLOTS}, IN THAT ORDER. 
Rscript $SCRIPT_DIR/TG.3D_runPQLseq.r ${SGE_TASK_ID} ${NSLOTS} TP1 MS

# Unload module
module unload R/4.2.2
