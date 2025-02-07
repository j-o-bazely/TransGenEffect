####                    TG.2_GlobalMethAnalysis                             ####
## Created: December 5th 2024
## By : James Bazely
## Context: This script analyses shared CpG between the three key treatment groups

#### ENV PREP ####
## Libraries
library(dplyr)
library(methylKit)
library(UpSetR)
library(ComplexHeatmap)

## Set paths 
dataDir <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage'

## load functions
source('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_analysis/TG_scripts/TG.0_Functions.r')

#### ANALYSIS ####

###
# Questions to ask:
#   Compare CpG within treatments between time points
#   Compare between treatments within time point

####            Comparing CpG Within Treatments, between time points            ####

TimePointCompPlots <- function(group){
    ## Load Objects
    print('Loading objects...')
    TP1_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TP1/TG_TP1_%sMethObj75pc.RDS', group)))
    TP3_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TP3/TG_TP3_%sMethObj75pc.RDS', group)))
    TP6_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TP6/TG_TP6_%sMethObj75pc.RDS', group)))
    
    ## Create granges
    print('Creating Grange objects...')
    TP1_grange <- makeGrangeFromMethObj(TP1_obj)
    print('     TP1 done!')
    TP3_grange <- makeGrangeFromMethObj(TP3_obj)
    print('     TP3 done!')
    TP6_grange <- makeGrangeFromMethObj(TP6_obj)
    print('     TP6 done!')
    
    ## Create a GRanges object of sites present in all 3 timepoints
    unionGRange <- GenomicRanges::intersect(TP1_grange, TP3_grange) %>%
      GenomicRanges::intersect(., TP6_grange)
    
    ## recreate Granges
    TP1_unionObj <- selectByOverlap(TP1_obj, unionGRange)
    TP3_unionObj <- selectByOverlap(TP3_obj, unionGRange)
    TP6_unionObj <- selectByOverlap(TP6_obj, unionGRange)
    
    ##  Box plot of mean methylation
    print('Creating percent methylation objects...')
    TP1perc <- data.frame(percMethylation(TP1_unionObj, rowids = TRUE)) %>%
      rownames_to_column(var = 'site') %>%
      mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
    print('     TP1 done!')
    
    TP3perc <- data.frame(percMethylation(TP3_unionObj, rowids = TRUE)) %>%
      rownames_to_column(var = 'site') %>%
      mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
    print('     TP3 done!')
    
    TP6perc <- data.frame(percMethylation(TP6_unionObj, rowids = TRUE)) %>%
      rownames_to_column(var = 'site') %>%
      mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
    print('     TP6 done!')

    
    print('Generating box plots...')
    boxplot <- ggplot() +
      ## Hatchling Exclusive
      geom_boxplot(
        data =  TP1perc,
        aes(x = sprintf("TP1 %s", group), y = mean_methylation),
        fill = "blue"
      ) + 
      geom_boxplot(
        data = TP3perc,
        aes(x = sprintf('TP3 %s', group), y = mean_methylation),
        fill = 'purple'
      ) +
      geom_boxplot(
        data = TP6perc,
        aes(x = sprintf('TP6 %s', group), y = mean_methylation),
        fill = 'red'
      )+
      ggtitle(sprintf('Mean Methylation of all CpG in %ss across three time points', group))
    
    ## Upset plot of shared sites within treatment across time point
    print('Generating UpSet plot...')
    upsetList <- list(TP1 = TP1_grange,
                      TP3 = TP3_grange,
                      TP6 = TP6_grange)
    
    m1 = make_comb_mat(upsetList)
    
    return(list(boxplot = boxplot, upset_plot = UpSet(m1)))
}

motherPlots <- TimePointCompPlots('mother')
hatchPlots <- TimePointCompPlots('hatch')
deepPlots <- TimePointCompPlots('deep')
shallowPlots <- TimePointCompPlots('shallow')


mean(colMeans(TP1perc[, !names(TP1perc) %in% "site"], na.rm = TRUE))
mean(colMeans(TP3perc[, !names(TP1perc) %in% "site"], na.rm = TRUE))
mean(colMeans(TP6perc[, !names(TP1perc) %in% "site"], na.rm = TRUE))



