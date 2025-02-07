####                    TG.3B_GeneticRelatednessMatrix                      ####
## Created: December 12th 2024
## By : James Bazely
## Context: This script creates a genetic relatedness matrix for each of the three
##          intersections being compared with PQLseq
## IMPORTANT: Order of matrix needs to match order of samples given to PQLseq ##

#### Prep environment #######
## Load packages ##
library(tidyverse)

## Set directories ##
DIR_PQL <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.3_data"
DIR_META <- "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_metadata"

                      #### TIMEPOINT 1 #### 
TP = 1

#### Matrix 1: SD ####
## Load metadata ##
DSmeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))

## Subset for selected timepoint ##
DSmeta2 <- DSmeta[DSmeta$Relocation == TP,]

## Start with a matrix of all zeroes ##
DSgenRelMat <- matrix(c(0), nrow = length(DSmeta2$Sample), ncol = length(DSmeta2$Sample))

# Add row and column names
rownames(DSgenRelMat) <- DSmeta2$Sample
colnames(DSgenRelMat) <- DSmeta2$Sample

DSrow_names <- rownames(DSgenRelMat)
DScol_names <- colnames(DSgenRelMat)

# If mother is the same for row_name and col_name, that value needs to be 0.5      <- THIS NEEDS TO BE FIXED
for (row_name in DSrow_names) {
  for(col_name in DScol_names){
    if(DSmeta2$Mother[DSmeta2$Sample == col_name] == DSmeta2$Mother[DSmeta2$Sample == row_name]){
      DSgenRelMat[col_name, row_name] <- 0.5
    }
  }
}

# If row name and column name are the same (i.e. same sample ID), fill cell with 1
for (row_name in DSrow_names) {
  if (row_name %in% DScol_names) {
    DSgenRelMat[row_name, row_name] <- 1 
  }
}  

# Change order to match input data for PQLseq
DSgenRelMatSorted <- DSgenRelMat[DSrow_names, DSrow_names]

pheatmap(DSgenRelMatSorted)

#### Matrix 2: MD ####
## Load metadata ##
MD_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MD_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MD_hatchMeta2 <- MD_hatchMeta[MD_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Deep')
MD_motherMeta2 <- MD_motherMeta[MD_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MDgenRelMat <- matrix(c(0), nrow = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID), ncol = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MD_hatchMeta2$Sample, MD_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MDgenRelMat) <- sampleVector
colnames(MDgenRelMat) <- sampleVector

MDrow_names <- rownames(MDgenRelMat)
MDcol_names <- colnames(MDgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MD_hatchMeta2$Sample
motherSampleVector <- MD_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == j] == MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i]){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j,i] = 0.5
    }
  }
}



## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i] == j){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j, i] = 0.5
    }
  }
}




## Check if rowname matches colname
for(i in rownames(MDgenRelMat)){
  for(j in colnames(MDgenRelMat)){
    if (i == j){
      MDgenRelMat[i, j] = 1
    }
  }
}

## Sort samples to match pqlSEQ order ##
MDgenRelMatSorted <- MDgenRelMat[MDrow_names, MDrow_names]

pheatmap(MDgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
         )


#### Matrix 3: MS ####
## Load metadata ##
MS_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MS_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MS_hatchMeta2 <- MS_hatchMeta[MS_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Shallow')
MS_motherMeta2 <- MS_motherMeta[MS_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MSgenRelMat <- matrix(c(0), nrow = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID), ncol = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MS_hatchMeta2$Sample, MS_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MSgenRelMat) <- sampleVector
colnames(MSgenRelMat) <- sampleVector

MSrow_names <- rownames(MSgenRelMat)
MScol_names <- colnames(MSgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MS_hatchMeta2$Sample
motherSampleVector <- MS_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == j] == MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i]){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j,i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)

## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i] == j){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j, i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)


## Check if rowname matches colname
for(i in rownames(MSgenRelMat)){
  for(j in colnames(MSgenRelMat)){
    if (i == j){
      MSgenRelMat[i, j] = 1
    }
  }
}

heatmap(MSgenRelMat)

## Sort samples to match pqlSEQ order ##
MSgenRelMatSorted <- MSgenRelMat[MSrow_names, MSrow_names]

pheatmap(MSgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
)






#### Save TP1 objects ####
write.csv(as.data.frame(DSgenRelMatSorted), file.path(DIR_PQL, 'TP1/TP1_DS_relatednessMatrix.csv'))
write.csv(as.data.frame(MDgenRelMatSorted), file.path(DIR_PQL, 'TP1/TP1_MD_relatednessMatrix.csv'))
write.csv(as.data.frame(MSgenRelMatSorted), file.path(DIR_PQL, 'TP1/TP1_MS_relatednessMatrix.csv'))

                          #### TIMEPOINT 3 ####
TP = 3

#### Matrix 1: SD ####
## Load metadata ##
DSmeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))

## Subset for selected timepoint ##
DSmeta2 <- DSmeta[DSmeta$Relocation == TP,]

## Start with a matrix of all zeroes ##
DSgenRelMat <- matrix(c(0), nrow = length(DSmeta2$Sample), ncol = length(DSmeta2$Sample))

# Add row and column names
rownames(DSgenRelMat) <- DSmeta2$Sample
colnames(DSgenRelMat) <- DSmeta2$Sample

DSrow_names <- rownames(DSgenRelMat)
DScol_names <- colnames(DSgenRelMat)

# If mother is the same for row_name and col_name, that value needs to be 0.5      <- THIS NEEDS TO BE FIXED
for (row_name in DSrow_names) {
  for(col_name in DScol_names){
    if(DSmeta2$Mother[DSmeta2$Sample == col_name] == DSmeta2$Mother[DSmeta2$Sample == row_name]){
      DSgenRelMat[col_name, row_name] <- 0.5
    }
  }
}

# If row name and column name are the same (i.e. same sample ID), fill cell with 1
for (row_name in DSrow_names) {
  if (row_name %in% DScol_names) {
    DSgenRelMat[row_name, row_name] <- 1 
  }
}  

# Change order to match input data for PQLseq
DSgenRelMatSorted <- DSgenRelMat[DSrow_names, DSrow_names]

pheatmap(DSgenRelMatSorted)



#### Matrix 2: MD ####
## Load metadata ##
MD_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MD_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MD_hatchMeta2 <- MD_hatchMeta[MD_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Deep')
MD_motherMeta2 <- MD_motherMeta[MD_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MDgenRelMat <- matrix(c(0), nrow = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID), ncol = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MD_hatchMeta2$Sample, MD_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MDgenRelMat) <- sampleVector
colnames(MDgenRelMat) <- sampleVector

MDrow_names <- rownames(MDgenRelMat)
MDcol_names <- colnames(MDgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MD_hatchMeta2$Sample
motherSampleVector <- MD_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == j] == MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i]){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j,i] = 0.5
    }
  }
}

heatmap(MDgenRelMat)

## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i] == j){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j, i] = 0.5
    }
  }
}

heatmap(MDgenRelMat)


## Check if rowname matches colname
for(i in rownames(MDgenRelMat)){
  for(j in colnames(MDgenRelMat)){
    if (i == j){
      MDgenRelMat[i, j] = 1
    }
  }
}

heatmap(MDgenRelMat)

## Sort samples to match pqlSEQ order ##
MDgenRelMatSorted <- MDgenRelMat[MDrow_names, MDrow_names]

pheatmap(MDgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
)


#### Matrix 3: MS ####
## Load metadata ##
MS_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MS_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MS_hatchMeta2 <- MS_hatchMeta[MS_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Shallow')
MS_motherMeta2 <- MS_motherMeta[MS_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MSgenRelMat <- matrix(c(0), nrow = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID), ncol = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MS_hatchMeta2$Sample, MS_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MSgenRelMat) <- sampleVector
colnames(MSgenRelMat) <- sampleVector

MSrow_names <- rownames(MSgenRelMat)
MScol_names <- colnames(MSgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MS_hatchMeta2$Sample
motherSampleVector <- MS_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == j] == MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i]){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j,i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)

## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i] == j){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j, i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)


## Check if rowname matches colname
for(i in rownames(MSgenRelMat)){
  for(j in colnames(MSgenRelMat)){
    if (i == j){
      MSgenRelMat[i, j] = 1
    }
  }
}

heatmap(MSgenRelMat)

## Sort samples to match pqlSEQ order ##
MSgenRelMatSorted <- MSgenRelMat[MSrow_names, MSrow_names]

pheatmap(MSgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
)

#### Save TP3 objects ####
write.csv(as.data.frame(DSgenRelMatSorted), file.path(DIR_PQL, 'TP3/TP3_DS_relatednessMatrix.csv'))
write.csv(as.data.frame(MDgenRelMatSorted), file.path(DIR_PQL, 'TP3/TP3_MD_relatednessMatrix.csv'))
write.csv(as.data.frame(MSgenRelMatSorted), file.path(DIR_PQL, 'TP3/TP3_MS_relatednessMatrix.csv'))



                              #### TIMEPOINT 6 ####
TP = 6

#### Matrix 1: SD ####
## Load metadata ##
DSmeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))

## Subset for selected timepoint ##
DSmeta2 <- DSmeta[DSmeta$Relocation == TP,]

## Start with a matrix of all zeroes ##
DSgenRelMat <- matrix(c(0), nrow = length(DSmeta2$Sample), ncol = length(DSmeta2$Sample))

# Add row and column names
rownames(DSgenRelMat) <- DSmeta2$Sample
colnames(DSgenRelMat) <- DSmeta2$Sample

DSrow_names <- rownames(DSgenRelMat)
DScol_names <- colnames(DSgenRelMat)

# If mother is the same for row_name and col_name, that value needs to be 0.5      <- THIS NEEDS TO BE FIXED
for (row_name in DSrow_names) {
  for(col_name in DScol_names){
    if(DSmeta2$Mother[DSmeta2$Sample == col_name] == DSmeta2$Mother[DSmeta2$Sample == row_name]){
      DSgenRelMat[col_name, row_name] <- 0.5
    }
  }
}

# If row name and column name are the same (i.e. same sample ID), fill cell with 1
for (row_name in DSrow_names) {
  if (row_name %in% DScol_names) {
    DSgenRelMat[row_name, row_name] <- 1 
  }
}  

# Change order to match input data for PQLseq
DSgenRelMatSorted <- DSgenRelMat[DSrow_names, DSrow_names]
heatmap(DSgenRelMatSorted, Rowv = NA, Colv = NA)
pheatmap(DSgenRelMatSorted)



#### Matrix 2: MD ####
## Load metadata ##
MD_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MD_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MD_hatchMeta2 <- MD_hatchMeta[MD_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Deep')
MD_motherMeta2 <- MD_motherMeta[MD_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MDgenRelMat <- matrix(c(0), nrow = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID), ncol = length(MD_hatchMeta2$Sample) + length(MD_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MD_hatchMeta2$Sample, MD_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MDgenRelMat) <- sampleVector
colnames(MDgenRelMat) <- sampleVector

MDrow_names <- rownames(MDgenRelMat)
MDcol_names <- colnames(MDgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MD_hatchMeta2$Sample
motherSampleVector <- MD_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == j] == MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i]){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j,i] = 0.5
    }
  }
}

heatmap(MDgenRelMat)

## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MD_hatchMeta2$Mother[MD_hatchMeta2$Sample == i] == j){
      MDgenRelMat[i, j] = 0.5
      MDgenRelMat[j, i] = 0.5
    }
  }
}

heatmap(MDgenRelMat)


## Check if rowname matches colname
for(i in rownames(MDgenRelMat)){
  for(j in colnames(MDgenRelMat)){
    if (i == j){
      MDgenRelMat[i, j] = 1
    }
  }
}

heatmap(MDgenRelMat)

## Sort samples to match pqlSEQ order ##
MDgenRelMatSorted <- MDgenRelMat[MDrow_names, MDrow_names]

pheatmap(MDgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
)


#### Matrix 3: MS ####
## Load metadata ##
MS_hatchMeta <- read.csv(file.path(DIR_META, 'Metadata_TSD_All_Hatchlings.csv'))
MS_motherMeta <- read.csv(file.path(DIR_META, 'TG_TPMotherMetadata.csv')) %>%
  dplyr::select('Female_ID', 'Time_Point') %>%
  distinct(Female_ID, .keep_all = TRUE)

## Subset for selected timepoint ##
MS_hatchMeta2 <- MS_hatchMeta[MS_hatchMeta$Relocation == TP,] %>%
  filter(Depth == 'Shallow')
MS_motherMeta2 <- MS_motherMeta[MS_motherMeta$Time_Point == TP,]

## Start with a matrix of all zeroes ##
MSgenRelMat <- matrix(c(0), nrow = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID), ncol = length(MS_hatchMeta2$Sample) + length(MS_motherMeta2$Female_ID))

## create rowname and colname vector
sampleVector <- c(MS_hatchMeta2$Sample, MS_motherMeta2$Female_ID)

## Add row and column names ##
rownames(MSgenRelMat) <- sampleVector
colnames(MSgenRelMat) <- sampleVector

MSrow_names <- rownames(MSgenRelMat)
MScol_names <- colnames(MSgenRelMat)


## Create separate list of hatchlings and mothers
hatchSampleVector <- MS_hatchMeta2$Sample
motherSampleVector <- MS_motherMeta2$Female_ID

## Check if hatchling samples have the same mother ##
for (i in hatchSampleVector) {
  print(i)
  for(j in hatchSampleVector){
    print(j)
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == j] == MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i]){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j,i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)

## Check if mother IS mother of hatchling ##
for(i in hatchSampleVector){
  for(j in motherSampleVector){
    if(MS_hatchMeta2$Mother[MS_hatchMeta2$Sample == i] == j){
      MSgenRelMat[i, j] = 0.5
      MSgenRelMat[j, i] = 0.5
    }
  }
}

heatmap(MSgenRelMat)


## Check if rowname matches colname
for(i in rownames(MSgenRelMat)){
  for(j in colnames(MSgenRelMat)){
    if (i == j){
      MSgenRelMat[i, j] = 1
    }
  }
}

heatmap(MSgenRelMat)

## Sort samples to match pqlSEQ order ##
MSgenRelMatSorted <- MSgenRelMat[MSrow_names, MSrow_names]

pheatmap(MSgenRelMatSorted,
         cluster_rows = FALSE, # Disable row clustering
         cluster_cols = FALSE, # Disable column clustering
)










#### Save TP6 objects ####
write.csv(as.data.frame(DSgenRelMatSorted), file.path(DIR_PQL, 'TP6/TP6_DS_relatednessMatrix.csv'))
write.csv(as.data.frame(MDgenRelMatSorted), file.path(DIR_PQL, 'TP6/TP6_MD_relatednessMatrix.csv'))
write.csv(as.data.frame(MSgenRelMatSorted), file.path(DIR_PQL, 'TP6/TP6_MS_relatednessMatrix.csv'))





# test #
testDF <- as.data.frame(DSgenRelMatSorted)
