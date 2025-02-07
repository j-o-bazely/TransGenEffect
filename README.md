# TransGenEffect


**Description**
- This is the code for my project looking at transgenerational effects between a mother and her hatchlings, how this relationship changes depending on incubation depth of hatchlings, and how this relationship changes across a warming nesting season. 



# Script Order and Contents #
**TG.0_functions.r**
-  To minimise code duplication, this pipeline takes a function-based approach where object preparation is performed by custom functions which work across timepoints and intersections.
-  All custom functions can be found in this script, grouped by the script they are used in. 
-  The start of every other script should source() to the location of TG.0_functions.r at the top to ensure functions are available for use.

**TG.1_objectPreparation.r**
- This script prepares 75% coverage objects for each key analysis group; mothers, hatchlings, deep hatchlings and shallow hatchlings.
- It also prepares these objects split into the three time points used within the split-clutch experiment, representing three different temperatures across the nesting season.
- Next, it creates the necessary GRanges objects to interrogate each intersection at each of the three time points.

**TG.2_GlobalMethAnalysis.r**
- This script analyses shared CpG between the three key treatment groups on the global methylation level



**-  -  -  -  -  -  -  -  -  -   TG.3_PQLseqAnalysis  -  -  -  -  -  -  -  -  -  - **
  - This section is split into multiple scripts

**TG.3A_prepMethDiffObj.r**
- This script generates diffmeth objects between each intersection, within each timepoint

**TG.3B_GeneticRelatednessMatrix.r**
- This script creates a genetic relatedness matrix for each of the three intersections being compared with PQLseq

**TG.3C_PrepPQLseqInputData.r**
- Prepares the methylKit uniteCov object for PQLseq's requirements                                 

**TG.3D_runPQLseq.r**
- Script for submission via Apocrita HPC. Runs PQLseq per chromosome

**TG.3E_runPQLseq_script.sh**
- Shell script for running TG.3D_runPQLseq.r

**TG.3F_mergeChrms_runSLMIM.r**
- This script runs the SLIM multiple testing algorithm, merges chromosomes into one object and combines with diffMeth object for downstream analysis

-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 




