## TG.0_functions ##
#### prepMethObj100pc ####
prepMethObj100pc <- function(metadata, group = c('mother', 'hatchling', 'deep', 'shallow'), tp = NA, nslots = nslots){
  
  # Start timer
  start.time <- Sys.time()
  
  # Path to sample methylation data
  SamplePathVector <- c()
  sampleIdList <- list()
  treatmentVector <- c()
  
  if(tp == 'all_three!!!'){
    if(group == 'mother'){
      uniqueFemale <- unique(metadata$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    
    if(group == 'hatchling'){
      hatchSamples <- unique(metadata$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- metadata %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- metadata %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
  }
  if(tp == '1'){
    tp1meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp1meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp1meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp1meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  if(tp == '3'){
    tp3meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp3meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp3meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp3meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  if(tp == '6'){
    tp6meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp6meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp6meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp6meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  # Running MethRead()
  methObj <- methRead(location = SamplePathList,
                      sample.id = sampleIdList,
                      treatment = treatmentVector,
                      context = 'CpG',
                      assembly = 'Yen2022',
                      pipeline="bismarkCoverage"
  )
  
  ## Coverage & Global Stats ##
  getMethylationStats(methObj[[1]], plot=TRUE)
  
  getCoverageStats(methObj[[1]], plot=TRUE)
  
  
  #### Data Cleaning & Object Merging                                           ####
  
  # Filter object by coverage (low coverage filter and >99.9th percentile in each sample)
  # Filter by 5X
  MethFilter5Cov <-filterByCoverage(methObj,
                                    lo.count=5,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)
  
  # Normalise coverage stats
  MethNormFilter5Cov=normalizeCoverage(MethFilter5Cov)
  
  # Merge samples together into a single united coverage object
  uniteMeth100pc <- methylKit::unite(MethNormFilter5Cov, destrand = FALSE, mc.cores = nslots)
  
  # Return this as the final object
  return(uniteMeth100pc)
}
#### prepMethObj75pc ####
prepMethObj75pc <- function(metadata, group = c('mother', 'hatchling', 'deep', 'shallow'), tp = NA, nslots = nslots){
  
  # Start timer
  start.time <- Sys.time()
  
  
  
  # Path to sample methylation data
  SamplePathVector <- c()
  sampleIdList <- list()
  treatmentVector <- c()
  
  if(tp == 'all_three!!!'){
    if(group == 'mother'){
      uniqueFemale <- unique(metadata$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    
    if(group == 'hatchling'){
      hatchSamples <- unique(metadata$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- metadata %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- metadata %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
  }
  if(tp == '1'){
    tp1meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp1meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp1meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp1meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  if(tp == '3'){
    tp3meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp3meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp3meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp3meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  if(tp == '6'){
    tp6meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'mother'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      for(i in 1:length(uniqueFemale)){
        treatmentVector <- append(treatmentVector, 1) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- uniqueFemale[i]
      }
    }
    if(group == 'hatchling'){
      hatchSamples <- unique(tp6meta$Sample)
      for(i in 1:length(hatchSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          hatchSamples
        )
      }
      for(i in 1:length(hatchSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- hatchSamples[i]
      }
    }
    if(group == 'deep'){
      deepSamples <- tp6meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      for(i in 1:length(deepSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- deepSamples[i]
        
      }
    }
    if(group == 'shallow'){
      shallowSamples <- tp6meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples
        )
      }
      for(i in 1:length(shallowSamples)){
        treatmentVector <- append(treatmentVector, 2) # 1 = mother | 2 = hatchling
        sampleIdList[i] <- shallowSamples[i]
      }
    }
    
    
    ## List of sample IDs
    SamplePathList <- list()
    
    for (i in 1:length(SamplePathVector)) {
      SamplePathList[i] <- SamplePathVector[i]
    }
    
  }
  
  # Running MethRead()
  methObj <- methRead(location = SamplePathList,
                      sample.id = sampleIdList,
                      treatment = treatmentVector,
                      context = 'CpG',
                      assembly = 'Yen2022',
                      pipeline="bismarkCoverage"
  )
  
  ## Coverage & Global Stats ##
  #getMethylationStats(methObj[[1]], plot=TRUE)
  
  #getCoverageStats(methObj[[1]], plot=TRUE)
  
  
  #### Data Cleaning & Object Merging                                           ####
  
  # Filter object by coverage (low coverage filter and >99.9th percentile in each sample)
  # Filter by 5X
  MethFilter5Cov <-filterByCoverage(methObj,
                                    lo.count=5,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)
  
  # Normalise coverage stats
  MethNormFilter5Cov=normalizeCoverage(MethFilter5Cov)
  
  # Calculate 75% of samples for min.per.group coverage
  minSamples <- as.integer(ceiling(length(sampleIdList)*0.75)) ## Gone conservative for now
  
  # Merge samples together into a single united coverage object
  uniteMeth75pc <- methylKit::unite(MethNormFilter5Cov, destrand = FALSE, min.per.group = minSamples, mc.cores = nslots)
  
  # Return this as the final object
  return(uniteMeth75pc)
}

#### makeGrangeFromMethObj ####

makeGrangeFromMethObj <- function(methObj){
  percObj <- percMethylation(methObj)
  methObj$original_pos <- methObj$start
  methObj$original_chrom <- methObj$chr
  
  objPos <- paste(methObj$chr, methObj$end, methObj$original_pos, methObj$original_chrom)
  
  objDF <- data.frame(chr=sapply(strsplit(objPos, " "), `[`, 1),
                      start=sapply(strsplit(objPos, " "), `[`, 2),
                      end=sapply(strsplit(objPos, " "), `[`, 2),
                      original_pos=sapply(strsplit(objPos, " "), `[`, 3),
                      original_chrom=sapply(strsplit(objPos, " "), `[`, 4))
  
  objGrange <- makeGRangesFromDataFrame(objDF, keep.extra.columns = TRUE)
  return(objGrange)
}





#### PrepGlobMethObj100pc ####
prepGlobMethObj100pc <- function(tp){
  mother  <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_motherMethObj100pc.RDS', tp, tp)))
  hatch   <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_hatchMethObj100pc.RDS', tp, tp)))
  deep    <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_deepMethObj100pc.RDS', tp, tp)))
  shallow <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_shallowMethObj100pc.RDS', tp, tp)))
  
  #### ANALYSIS ####
  ####    GRANGES ###                                                             ####
  Mgrange <- makeGrangeFromMethObj(mother)
  Hgrange <- makeGrangeFromMethObj(hatch)
  Dgrange <- makeGrangeFromMethObj(deep)
  Sgrange <- makeGrangeFromMethObj(shallow)
  
  ####      OVERLAP GRANGES                                                       ####
  MH <- subsetByOverlaps(Mgrange, Hgrange)
  MS <- subsetByOverlaps(Mgrange, Sgrange)
  MD <- subsetByOverlaps(Mgrange, Dgrange)
  SD <- subsetByOverlaps(Sgrange, Dgrange)
  SDM <- subsetByOverlaps(SD, Mgrange)
  
  ####        SAVE OBJECTS                                                        ####
  saveRDS(MH, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_MHgrange.RDS', tp, tp)))
  saveRDS(MS, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_MSgrange.RDS', tp, tp)))
  saveRDS(MD, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_MDgrange.RDS', tp, tp)))
  saveRDS(SD, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_SDgrange.RDS', tp, tp)))
  saveRDS(SDM, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_SDMgrange.RDS', tp, tp)))

  ####          EXCLUSIVE INTERSECTIONS ####       ## Note for me: Don't overcomplicate, just create exclusive intersections for now, setdiff can easily be run in downstream analysis to create objects needed
  MS_exclusive <- setdiff(MS, SDM)
  SD_exclusive <- setdiff(SD, SDM)
  MD_exclusive <- setdiff(MD, SDM)
  
  ####            SAVE OBJECTS                                                        ####
  saveRDS(MS_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_MSgrange_exclusive.RDS', tp, tp)))
  saveRDS(SD_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_SDgrange_exclusive.RDS', tp, tp)))
  saveRDS(MD_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_%s/TG_%s_MDgrange_exclusive.RDS', tp, tp)))
  }

#### TimePointCompPlots ####
TimePointCompPlots <- function(group){
  ## Load Objects
  print('Loading objects...')
  TP1_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP1_%sMethObj100pc.RDS', group)))
  TP3_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP3/TG_TP3_%sMethObj100pc.RDS', group)))
  TP6_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP6/TG_TP6_%sMethObj100pc.RDS', group)))
  
  print('Creating Grange objects...')
  TP1_grange <- makeGrangeFromMethObj(TP1_obj)
  print('     TP1 done!')
  TP3_grange <- makeGrangeFromMethObj(TP3_obj)
  print('     TP3 done!')
  TP6_grange <- makeGrangeFromMethObj(TP6_obj)
  print('     TP6 done!')
  
  ##  Box plot of mean methylation
  print('Creating percent methylation objects...')
  TP1perc <- data.frame(percMethylation(TP1_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP1 done!')
  
  TP3perc <- data.frame(percMethylation(TP3_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP3 done!')
  
  TP6perc <- data.frame(percMethylation(TP6_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP6 done!')
  
  ### NOTE: I realise this will hold more sway if the sites are identical: Make sure to test this after doing the upset plot
  print('Generating box plots...')
  boxplot <- ggplot() +
    ## Hatchling Exclusive
    geom_boxplot(
      data = TP1perc,
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
  
  ## Subset for shared sites and re-plot the box plot
  
  # Load SDM objects for time points
  TP1_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP1_SDMgrange', group)))
  TP3_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP3_SDMgrange', group)))
  TP6_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP6_SDMgrange', group)))
  
  # Subset methObjs by these gRanges
  TP1_obj_SDM <- selectByOverlap(TP1_SDMgrange, TP1_obj)
  TP3_obj_SDM <- selectByOverlap(TP3_SDMgrange, TP3_obj)
  TP6_obj_SDM <- selectByOverlap(TP6_SDMgrange, TP6_obj)
  
  # 
  

  return(list(boxplot = boxplot, upset_plot = UpSet(m1)))
  
  ## Subset for shared sites and re-plot the box plot for these
  
}


#### prepGlobalMethObj75pc ####
prepGlobMethObj75pc <- function(tp){
  mother  <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_motherMethObj75pc.RDS', tp, tp)))
  hatch   <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_hatchMethObj75pc.RDS', tp, tp)))
  deep    <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_deepMethObj75pc.RDS', tp, tp)))
  shallow <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_shallowMethObj75pc.RDS', tp, tp)))
  
  #### ANALYSIS ####
  ####    GRANGES ###                                                             ####
  Mgrange <- makeGrangeFromMethObj(mother)
  Hgrange <- makeGrangeFromMethObj(hatch)
  Dgrange <- makeGrangeFromMethObj(deep)
  Sgrange <- makeGrangeFromMethObj(shallow)
  
  ####      OVERLAP GRANGES                                                       ####
  MH <- subsetByOverlaps(Mgrange, Hgrange)
  MS <- subsetByOverlaps(Mgrange, Sgrange)
  MD <- subsetByOverlaps(Mgrange, Dgrange)
  DS <- subsetByOverlaps(Sgrange, Dgrange)
  SDM <- subsetByOverlaps(DS, Mgrange)
  
  ####        SAVE OBJECTS                                                        ####
  saveRDS(MH, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_MHgrange75pc.RDS', tp, tp)))
  saveRDS(MS, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_MSgrange75pc.RDS', tp, tp)))
  saveRDS(MD, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_MDgrange75pc.RDS', tp, tp)))
  saveRDS(DS, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_DSgrange75pc.RDS', tp, tp)))
  saveRDS(SDM, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_SDMgrange75pc.RDS', tp, tp)))
  
  ####          EXCLUSIVE INTERSECTIONS ####       ## Note for me: Don't overcomplicate, just create exclusive intersections for now, setdiff can easily be run in downstream analysis to create objects needed
  MS_exclusive <- setdiff.Vector(MS, SDM)
  DS_exclusive <- setdiff.Vector(DS, SDM)
  MD_exclusive <- setdiff.Vector(MD, SDM)
  
  ####            SAVE OBJECTS                                                        ####
  saveRDS(MS_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_MSgrange_exclusive75pc.RDS', tp, tp)))
  saveRDS(DS_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_DSgrange_exclusive75pc.RDS', tp, tp)))
  saveRDS(MD_exclusive, file.path(dataDir, sprintf('TG_analysisData/TG.1_data/%s/TG_%s_MDgrange_exclusive75pc.RDS', tp, tp)))
}

#### TimePointCompPlots ####
TimePointCompPlots <- function(group){
  ## Load Objects
  print('Loading objects...')
  TP1_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP1_%sMethObj75pc.RDS', group)))
  TP3_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP3/TG_TP3_%sMethObj75pc.RDS', group)))
  TP6_obj <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP6/TG_TP6_%sMethObj75pc.RDS', group)))
  
  print('Creating Grange objects...')
  TP1_grange <- makeGrangeFromMethObj(TP1_obj)
  print('     TP1 done!')
  TP3_grange <- makeGrangeFromMethObj(TP3_obj)
  print('     TP3 done!')
  TP6_grange <- makeGrangeFromMethObj(TP6_obj)
  print('     TP6 done!')
  
  ##  Box plot of mean methylation
  print('Creating percent methylation objects...')
  TP1perc <- data.frame(percMethylation(TP1_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP1 done!')
  
  TP3perc <- data.frame(percMethylation(TP3_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP3 done!')
  
  TP6perc <- data.frame(percMethylation(TP6_obj, rowids = TRUE)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(mean_methylation = rowMeans(.[-1], na.rm = TRUE))
  print('     TP6 done!')
  
  ### NOTE: I realise this will hold more sway if the sites are identical: Make sure to test this after doing the upset plot
  print('Generating box plots...')
  boxplot <- ggplot() +
    ## Hatchling Exclusive
    geom_boxplot(
      data = TP1perc,
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
  
  ## Subset for shared sites and re-plot the box plot
  
  # Load SDM objects for time points
  TP1_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP1_SDMgrange', group)))
  TP3_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP3_SDMgrange', group)))
  TP6_SDMgrange <- readRDS(file.path(dataDir, sprintf('TG_analysisData/TG.1_data/TG_TP1/TG_TP6_SDMgrange', group)))
  
  # Subset methObjs by these gRanges
  TP1_obj_SDM <- selectByOverlap(TP1_SDMgrange, TP1_obj)
  TP3_obj_SDM <- selectByOverlap(TP3_SDMgrange, TP3_obj)
  TP6_obj_SDM <- selectByOverlap(TP6_SDMgrange, TP6_obj)
  
  # 
  
  
  return(list(boxplot = boxplot, upset_plot = UpSet(m1)))
  
  ## Subset for shared sites and re-plot the box plot for these
  
}

#### PrepIntersectObj100pc ####
PrepIntersectObj100pc <- function(metadata, group = c('DS', 'MD', 'MS'), tp = NA,  nslots = nslots){
  
  # Path to sample methylation data
  SamplePathVector <- c()
  sampleIdList <- list()
  treatmentVector <- c()
  treatmentVectorStr <- c()
  
  # Start timer
  start.time <- Sys.time()
  
  if(tp == '1'){
    tp1meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp1meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp1meta$Depth[tp1meta$Sample == samples[i]]) # 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp1meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp1meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples,
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }
  if(tp == '3'){
    tp3meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp3meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp3meta$Depth[tp3meta$Sample == samples[i]])# 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp3meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp3meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples,
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }
  if(tp == '6'){
    tp6meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp6meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp6meta$Depth[tp6meta$Sample == samples[i]])# 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp6meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp6meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples,
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }

  # Running MethRead()
  methObj <- methRead(location = SamplePathList,
                      sample.id = sampleIdList,
                      treatment = treatmentVector,
                      context = 'CpG',
                      assembly = 'Yen2022',
                      pipeline="bismarkCoverage"
  )
  
  ## Coverage & Global Stats ##
  #getMethylationStats(methObj[[1]], plot=TRUE)
  
  #getCoverageStats(methObj[[1]], plot=TRUE)
  
  
  #### Data Cleaning & Object Merging                                           ####
  
  # Filter object by coverage (low coverage filter and >99.9th percentile in each sample)
  # Filter by 5X
  MethFilter5Cov <-filterByCoverage(methObj,
                                    lo.count=5,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)
  
  # Normalise coverage stats
  MethNormFilter5Cov=normalizeCoverage(MethFilter5Cov)
  
  # Merge samples together into a single united coverage object
  uniteMeth100pc <- methylKit::unite(MethNormFilter5Cov, destrand = FALSE, mc.cores = nslots)

  # Finish timing
  finish.time <- Sys.time()
  print(finish.time - start.time)
  
  # Return this as the final object
  return(uniteMeth100pc)
}

#### PrepIntersectObj75pc ####
PrepIntersectObj75pc <- function(metadata, group = c('DS', 'MD', 'MS'), tp = NA,  nslots = nslots){
  
  # Path to sample methylation data
  SamplePathVector <- c()
  sampleIdList <- list()
  treatmentVector <- c()
  treatmentVectorStr <- c()
  
  # Start timer
  start.time <- Sys.time()
  
  if(tp == '1'){
    tp1meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp1meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp1meta$Depth[tp1meta$Sample == samples[i]])# 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
      
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp1meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp1meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp1meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples,
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }
  if(tp == '3'){
    tp3meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp3meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp3meta$Depth[tp3meta$Sample == samples[i]])# 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp3meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp3meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp3meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }
  if(tp == '6'){
    tp6meta <- metadata %>%
      filter(Relocation == tp)
    if(group == 'DS'){
      samples <- tp6meta %>%
        pull(Sample)
      for(i in 1:length(samples)){
        SamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          samples
        )
      }
      for(i in 1:length(samples)){
        treatmentVectorStr <- append(treatmentVectorStr, tp6meta$Depth[tp6meta$Sample == samples[i]])# 3 = deep  | 4 = shallow
        treatmentVector <- ifelse(treatmentVectorStr == 'Deep', 3, 4)
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MS'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      shallowSamples <- tp6meta %>%
        filter(Depth== 'Shallow') %>%
        pull(Sample)
      for(i in 1:length(shallowSamples)){
        shallowSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          shallowSamples,
          shallowSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, shallowSamples)
      SamplePathVector <- c(femaleSamplePathVector, shallowSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
    if(group == 'MD'){
      uniqueFemale <- unique(tp6meta$Mother)
      for(i in 1:length(uniqueFemale)){
        femaleSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_mothers/%s.CpG_merged.cov.gz",
          uniqueFemale
        )
      }
      deepSamples <- tp6meta %>%
        filter(Depth== 'Deep') %>%
        pull(Sample)
      for(i in 1:length(deepSamples)){
        deepSamplePathVector <- sprintf(
          "/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_hatchlings/%s.CpG_merged.cov.gz",
          deepSamples,
          deepSamples
        )
      }
      
      # Merge sample lists
      samples <- c(uniqueFemale, deepSamples)
      SamplePathVector <- c(femaleSamplePathVector, deepSamplePathVector)
      
      for(i in 1:length(samples)){
        treatmentVector <- append(treatmentVector, ifelse(grepl('SLL', samples[i]), 1, 4))
        sampleIdList[i] <- samples[i]
      }
      ## List of sample IDs
      SamplePathList <- list()
      
      for (i in 1:length(SamplePathVector)) {
        SamplePathList[i] <- SamplePathVector[i]
      }
    }
  }
  
  # Running MethRead()
  methObj <- methRead(location = SamplePathList,
                      sample.id = sampleIdList,
                      treatment = treatmentVector,
                      context = 'CpG',
                      assembly = 'Yen2022',
                      pipeline="bismarkCoverage"
  )
  
  ## Coverage & Global Stats ##
  #getMethylationStats(methObj[[1]], plot=TRUE)
  
  #getCoverageStats(methObj[[1]], plot=TRUE)
  
  
  #### Data Cleaning & Object Merging                                           ####
  
  # Filter object by coverage (low coverage filter and >99.9th percentile in each sample)
  # Filter by 5X
  MethFilter5Cov <-filterByCoverage(methObj,
                                    lo.count=5,
                                    lo.perc=NULL,
                                    hi.count=NULL,
                                    hi.perc=99.9)
  
  # Normalise coverage stats
  MethNormFilter5Cov=normalizeCoverage(MethFilter5Cov)
  
  # Calculate 75% of samples for min.per.group coverage
    # Create table of treatment Vector
  trtTable <- table(treatmentVector)
    # Select smallest result
  minGroupSize <- min(trtTable)
  minGroupSize
    # Find 75% of smallest result
  minSamples <- as.integer(ceiling(minGroupSize*0.75)) ## Gone conservative for now
  minSamples
  
  # Merge samples together into a single united coverage object
  uniteMeth75pc <- methylKit::unite(MethNormFilter5Cov, destrand = FALSE, min.per.group = minSamples, mc.cores = nslots)
  
  # Finish timing
  finish.time <- Sys.time()
  print(finish.time - start.time)
  
  # Return this as the final object
  return(uniteMeth75pc)
}

#### PrepNumCsFun ####
PrepNumCsFun <- function(list_of_dfs, uniteCov) {
  for (i in seq_along(list_of_dfs)) {
    # Delete chr column
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)
    
    # Remove columns containing "coverage" in their name
    col_names <- colnames(list_of_dfs[[i]]) # Delete all coverage columns
    cols_to_delete <- grep("^coverage", col_names)
    list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete] # Remove identified column
    
    # Rename coverage column names with just number, so it matches between the two input files
    # Could use sample.id, but easier to spot anything out of order with numbers
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% rename_with(~str_remove(., "^numCs"))
    #colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
    colnames(list_of_dfs[[i]]) <- paste(colnames(list_of_dfs[[i]]), uniteCov@sample.ids, sep = "_")
    
    
    
    
    # # Save row names for adding back in, as transformation step removes them (if performing)
    # row_names <- rownames(list_of_dfs[[i]])
    
    # # If doing count transformation: add +1 to the no. of methylated reads for every value
    # list_of_dfs[[i]] <- as.data.frame(lapply(list_of_dfs[[i]], function(x) if(is.numeric(x)) x + 1 else x))
    
    # # Add row names back in
    # rownames(list_of_dfs[[i]]) <- row_names
  }
  return(list_of_dfs)
}

#### PrepCoverageFun ####
PrepCoverageFun <- function(list_of_dfs, uniteCov) {
  for (i in seq_along(list_of_dfs)) {
    # Delete chr column
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% dplyr::select(-chr)
    
    # Delete all numCs columns
    col_names <- colnames(list_of_dfs[[i]]) # Delete all coverage columns
    cols_to_delete <- grep("^numCs", col_names)
    list_of_dfs[[i]] <- list_of_dfs[[i]][, -cols_to_delete] # Remove identified column
    
    # Rename coverage column names with just number, so it matches between the two input files
    # Could use sample.id, but easier to spot anything out of order with numbers
    list_of_dfs[[i]] <- list_of_dfs[[i]] %>% rename_with(~str_remove(., "coverage"))
    #colnames(list_of_dfs[[i]]) <- paste0("ind", colnames(list_of_dfs[[i]]))
    colnames(list_of_dfs[[i]]) <- paste(colnames(list_of_dfs[[i]]), uniteCov@sample.ids, sep = "_")
    
    # # Save row names for adding back in, as next step removes them
    #row_names <- rownames(list_of_dfs[[i]])
    
    # # If doing count transformation: add +2 to the no. of reads for every value
    # list_of_dfs[[i]] <- as.data.frame(lapply(list_of_dfs[[i]], function(x) if(is.numeric(x)) x + 2 else x))
    
    # # Add row names back in
    # rownames(list_of_dfs[[i]]) <- row_names
  }
  return(list_of_dfs)
}

#### PrepInputObjs ####
PrepInputObjs <- function(DIR_UniteCov, DIR_CHROM, DIR_PQL, INSECT, TP) {
  # Construct filenames based on INSECT and TP arguments
  uniteCov_file <- file.path(DIR_UniteCov, paste0("", TP, "/TG_", TP, "_", INSECT, "methObj75pc.RDS"))
  
  ## Load uniteCov object
  uniteCov <- readRDS(uniteCov_file)
  
  ## Pull out treatment vector
  trtVec <- uniteCov@treatment
  
  ## Convert uniteCov to dataframe
  uniteCov_df <- as.data.frame(getData(uniteCov))
  
  ## Rename rows to site names
  rownames(uniteCov_df) <- paste(uniteCov_df$chr, uniteCov_df$start, sep = ".")
  
  ## Remove site, start, end and strand columns (keep chr column)
  uniteCov_df <- uniteCov_df %>% dplyr::select(-start, -end, -strand)
  
  ## Delete ALL numTs columns
  col_names <- colnames(uniteCov_df)
  cols_to_delete <- grep("^numTs", col_names)
  uniteCov_df <- uniteCov_df[, -cols_to_delete]
  
  ## Split into separate DFs per chromosome
  uniteCov_df$chr <- as.factor(uniteCov_df$chr)
  list_chr_df <- split(uniteCov_df, uniteCov_df$chr)
  
  ## Prep RawCountDataSet input file
  list_chr_df_numCs <- PrepNumCsFun(list_chr_df, uniteCov)
  
  ## Prep LibSize input file
  list_chr_df_coverage <- PrepCoverageFun(list_chr_df, uniteCov)
  
  ## Save each chromosome to separate RDS for numCs and coverage
  lapply(names(list_chr_df_numCs), function(i) {
    saveRDS(list_chr_df_numCs[[i]], file.path(DIR_CHROM, paste0(TP, "/", TP , '_', INSECT, '_ChromObjs/', TP , '_', INSECT, "_Input_NumCs_", i, ".RDS")))
  })
  
  lapply(names(list_chr_df_coverage), function(i) {
    saveRDS(list_chr_df_coverage[[i]], file.path(DIR_CHROM, paste0(TP, "/", TP , '_', INSECT, '_ChromObjs/', TP , '_', INSECT, "_Input_Coverage_", i, ".RDS")))
  })
  
  ## Save full list objects
  saveRDS(list_chr_df_coverage, file.path(DIR_PQL, paste0(TP, "/", TP , '_', INSECT, "_list_chr_df_coverage.RDS")))
  saveRDS(list_chr_df_numCs, file.path(DIR_PQL, paste0(TP, "/", TP , '_', INSECT, "_list_chr_df_numCs.RDS")))
  
  ## Save treatment vector
  writeLines(as.character(trtVec), file.path(DIR_PQL, paste0(TP, "/", TP , '_', INSECT, "_treatmentVector.txt")))
  
  return(list(numCs = list_chr_df_numCs, coverage = list_chr_df_coverage, treatment = trtVec))
}

### TG.3A_prepMethDiffObj.r ###
prepMethDiffObj <- function(TP, INTER){
  ## Load object
  intersectObj <- readRDS(file.path(paste(DIR_DATA, '/TG_analysisData/TG.1_data/', TP, '/TG_', TP, '_', INTER, 'methObj75pc.RDS', sep = '')))
  
  ## Generate methDiff object
  diffMethObj <- calculateDiffMeth(intersectObj, mc.cores = 10)
  
  ## Save methDiff object
  saveRDS(diffMethObj, file.path(paste(DIR_DATA, '/TG_analysisData/TG.3_data/', TP, '/TG_', TP, '_', INTER, 'diffMethObj75pc.RDS', sep = '')))
}



#### SLIMfunc ####
#### F1 ####
f1 <- function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)*100}
SLIMfunc<-function(rawp,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)
{
  ####################
  m <- length(rawp) 
  
  ########################
  alpha_mtx=NULL;#
  pi0s_est_COM=NULL;
  SzCluster_mtx=NULL;
  P_pi1_mtx=NULL;
  pi1_act_mtx=NULL;
  Dist_group=NULL;
  Num_mtx=NULL;
  Gamma=NULL;
  PI0=NULL;
  
  #############
  ##observed points
  lambda_ga=seq(0,1,0.001);
  gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
  Gamma=c(Gamma,gamma_ga);
  alpha_mtx=c(alpha_mtx,gamma_ga[which(lambda_ga==0.05)]);
  
  
  ###esimation
  pi0_mtx=NULL;
  x.axis=NULL;
  y.axis=NULL;
  itv=(1-STA)/Divi;
  for (i in 1:Divi)##10 for uniform data
  {
    cutoff=STA+(i/Divi)*(1-STA);##10 for uniform data
    lambda=seq(cutoff-itv,cutoff,itv/10);
    gamma_mtx=sapply(lambda,f1,rawp=rawp);
    LModel=lm(gamma_mtx~lambda);
    pi0_mtx=c(pi0_mtx,coefficients(LModel)[2]);
  }
  
  
  ##################################
  ########searching
  N_COM=NULL;
  N_rawp=NULL;
  maxFDR_mtx=NULL;
  quapoint_mtx=NULL; 
  if (B<=1) B=100;
  quapoint_mtx=seq(0.01,0.99,1/B);
  for (k in 1:length(quapoint_mtx))
  {
    qua_point=quapoint_mtx[k];
    
    pi0_combLR=min(quantile(pi0_mtx,qua_point),1);#mean(pi0_mtx);#median();# qua_point=0.78 for desreasing distribution;
    ##0.4 for uniform or normal distribution;
    pi0_est=pi0_combLR;
    
    ###########Calculate independent index of raw p vlaues
    PI0=rbind(PI0,pi0_mtx);
    
    pi0s_est_COM=c(pi0s_est_COM,pi0_est);
    ##Condition1
    P_pi1=sort(rawp)[max(length(rawp)*(1-pi0_est),1)];##
    P_pi1_mtx=c(P_pi1_mtx,P_pi1);
    
    pi0=pi0_est;
    if (is.null(Pz)) Pz=0.05;
    maxFDR=Pz*pi0/(1-(1-Pz)*pi0);
    maxFDR_mtx=c(maxFDR_mtx,maxFDR);
    
    qvalues_combLR=QValuesfun(rawp,pi0);
    qvalue_cf=maxFDR;
    selected=which(qvalues_combLR<qvalue_cf);
    Sel_qvalues_combLR=selected;
    
    pi1_act_mtx=c(pi1_act_mtx,length(Sel_qvalues_combLR)/length(rawp));
    N_COM=c(N_COM,list(Sel_qvalues_combLR));
    Num_mtx=c(Num_mtx,length(Sel_qvalues_combLR));
  }
  length(N_COM)
  length(quapoint_mtx)
  
  ####doing judging
  ##by max FDR
  pi1s_est_COM=1-pi0s_est_COM;
  Diff=sum(rawp<=Pz)/length(rawp)-pi1_act_mtx;
  
  ###
  loc=which.min(abs(Diff));
  Diff.loc=Diff[loc];
  selQuantile=quapoint_mtx[loc];
  pi0_Est=min(1,pi0s_est_COM[loc]);
  maxFDR.Pz=Pz*pi0_Est/(1-(1-Pz)*pi0_Est);
  
  if(Bplot)
  {
    windows();
    par(mfrow=c(1,2));
    hist(rawp,main="Histogram of p-value");
    gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
    plot(lambda_ga,gamma_ga,type="l",main="Relationship of p- and q value",xlab=expression(lambda),ylab=expression(gamma),cex.lab=1.45,cex.axis=1.42)
    #par(xaxp=c(0,1,10));
    #axis(1);
    ##qvalues
    qValues=QValuesfun(rawp,pi0=pi0_Est);
    gammaq_ga=sapply(lambda_ga,f1,rawp=qValues);
    lines(lambda_ga,gammaq_ga,col="blue",lwd=2,lty="dashed")
    abline(v=Pz,col="black",lwd=2,lty="dotdash")
    abline(v=maxFDR.Pz,col="blue",lwd=2,lty="dotdash")
    text(0.75,0.6,labels=paste("L=",round(abs(Diff.loc),4),sep=""));
    leg=list(bquote("CPD of p-value"),bquote("CPD of q-value"),bquote("Pmax"==.(Pz)),bquote("FDRmax"==.(round(maxFDR.Pz,2))));
    legend("bottomright",legend=as.expression(leg),lwd=2,lty=c("solid","dashed","dotdash","dotdash"),col=c("black","blue","black","blue"));
  }
  
  return(list(pi0_Est=pi0_Est,selQuantile=selQuantile));
}

#### QValuesFun ####
QValuesfun<-function(rawp,pi0)
{
  order_rawp=sort(rawp);
  qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
  temp=cummin(qvalues[seq(length(qvalues),1,-1)])
  qvalues=temp[seq(length(temp),1,-1)];
  qvalues=qvalues[order(order(rawp))]
}



