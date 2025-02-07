#### TG.3.F_mergeChrms_runSLIM.r ####   
## Created: February 3rd 2025                                                   ##
## By : James Bazely                                                            ##  
## Context: This script runs the SLIM multiple testing algorithm,               ##
## merges chromosomes into one object                                           ##
## and combines with diffMeth object for downstream analysis                    ##



## Load packages ##
library(tidyverse)
library(methylKit)
library(limma)


## Set paths ##
WORKDIR <- '/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/'
DATADIR <- paste0(WORKDIR, 'TG_dataStorage')
PQLDIR <- paste0(WORKDIR, 'TG_analysis/TG_scripts/TG.3_PQLseq')
OUTDIR <- paste0(DATADIR, '/TG_analysisData/TG.4_data/', TP)

### ###
### THIS WILL NEED TO BE FUNCTIONISED ###
### INITIAL WRITING AND TESTING WILL BE WITH TP1DS ###
### ###

preparePQLseqOutput <- function(TP, INT){
  
  # Set OUTDIR
  OUTDIR <- paste0(DATADIR, '/TG_analysisData/TG.4_data/', TP)
  
  print('#### MERGE CHROMOSOMES ####')
  list_fit <- list()
  list_dfs_fit <- list()
  
  ## loop through each file, adding to list of dfs
  ## Create list of chromosomes ##
  chrom_number <- 1:28
  chrom_number <- append(chrom_number, '0', after = 28)
  chroms <- c(paste0('chr', chrom_number))
  
  ## Load chromosome code list ##
  chrList <- read.csv(file.path(DATADIR, 'TG_metadata/chrom_names.txt'), header = FALSE)
  chrList$V1 <- gsub('SLK063_ragtag_', '', chrList$V1)
  
  ## Set file prefix
  prefix <- paste0(TP, '_', INT, '_', 'fit_')
  
  ## set CHROMDIR ##
  CHROMDIR <- paste0('/data/SBCS-EizaguirreLab/James_B/cleanPHD/transGen/TG_dataStorage/TG_analysisData/TG.3_data/', TP, '/', TP, '_', INT, '_Fit')  
  
  ## Parse chromosome list
  for (chr in chrList$V1) {
    print(chr)
    list_fit[[ chr ]] <- readRDS(file.path(CHROMDIR, paste0(prefix, chr, ".RDS")))
    list_dfs_fit[[ chr ]] <- list_fit[[ chr ]]
  }
  
  # Bind together into a data frame
  df_AllChrms <- do.call(rbind, list_dfs_fit)
  
  
  # % of NaNs for h2
  length(df_AllChrms$h2[grep("NaN", df_AllChrms$h2)])/nrow(df_AllChrms)*100 # 14.9% NAs <- Might be a produce of using 75% objects
  
  
  # For some reason, chr<n> has been added to start of every row name again -> remove
  rownames(df_AllChrms) <- gsub("^.*?\\.", "", rownames(df_AllChrms)) # Remove everything before first period
  
  
  #### FUNCTIONS: TO BE MOVED TO FUNCTIONS SHEET ####
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
  
  f1 <- function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)*100};
  
  ###
  # Define function to calculate q-values
  QValuesfun<-function(rawp,pi0)
  {
    order_rawp=sort(rawp);
    qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
    temp=cummin(qvalues[seq(length(qvalues),1,-1)])
    qvalues=temp[seq(length(temp),1,-1)];
    qvalues=qvalues[order(order(rawp))]
  }
  
  
  pvals <- df_AllChrms$pvalue
  length(pvals)
  
  # Run SLIM
  qvals = {QValuesfun(pvals,
                      SLIMfunc(pvals,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)$pi0_Est)
  }
  
  slim_output <- df_AllChrms
  slim_output$qvalue <- qvals
  
  
  ###### Checks ######
  
  ## How many failed to converge? ##
  nrow(slim_output[slim_output$converged == "FALSE",]) # 6,168 sites
  nrow(slim_output[slim_output$qvalue == 0,]) # 1,309 sites
  nrow(slim_output[slim_output$qvalue == 0 & slim_output$converged == "FALSE",]) # 1,309 sites
  # -> All sites that are q=0 also failed to converge -> represent algorithm failure and should be removed
  
  ## How many have NaN in sig rows? ##
  q0.05 <- slim_output[slim_output$qvalue < 0.05,] # 5,097 sites
  nrow(q0.05[q0.05$h2 == "NaN",]) # 110 sites
  # How many failed to converge in sig rows?
  nrow(q0.05[q0.05$converged == "FALSE",]) # 4747 sites
  # How many have q=0 in sig rows
  nrow(q0.05[q0.05$qvalue == 0,]) # 1309 sites
  nrow(q0.05[q0.05$qvalue == 0 & q0.05$converged == "FALSE",]) #  7,203 sites
  
  # Summary ##
  #slim_output <- slim_output[slim_output$converged == "TRUE",] # 1309 sites left
  summary(slim_output$qvalue)
  # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  # 0.0000817 0.9984423 0.9984423 0.9919617 0.9984423 0.9984423 
  
  ## Histogram ##
  hist(slim_output$pvalue)
  hist(slim_output$qvalue)
  
  
  print('#### Combine with diffMeth object ####')
  methdiff <- readRDS(paste0(DATADIR, '/TG_analysisData/TG.3_data/', TP, '/', TP, '_', INT, '_diffMethObj75pc.RDS')) # This has been saved as a large character? No idea.
  methdiff <- readRDS('/data/scratch/bty500/TP1DStestMethDiff.RDS')
  
  ## Format methylKit object ##
  methdiff_df <- getData(methdiff) %>%
    mutate('pos' = paste0(chr, '.', start)) %>%
    dplyr::select(pos, meth.diff)
  
  ## Format PQLseq output ##
  slim_output <- slim_output %>%
    rownames_to_column('pos')
  
  ## Join the two DFs by 
  methdiff_df <- left_join(methdiff_df, slim_output, by = 'pos')
  
  print('#### Filtering ####')
  ## A: Remove sites that failed to converge ##
  A <- methdiff_df[methdiff_df$converged == TRUE,]
  
  ## B: Remove sites with q=0 ##
  B <- A[A$qvalue != 0,]
  
  ## C: Remove sites with h2=NaN ##
  C <- B[!is.na(B$h2),]
  
  ## D: Remove sites that don't meet the significance threshold ##
  qThresh = 0.05
  D <- C[C$qvalue < qThresh,]
  
  ## E: Remove sites that don't meet the minimum % difference threshold ##
  percThresh = 10
  E <- D[abs(D$meth.diff) > percThresh,]


  print('#### Create percent Methylation object ####')
  ## Load uniteCov object ##
  unitecov <- readRDS(file.path(paste0(DATADIR, '/TG_analysisData/TG.1_data/', TP, '/TG_', TP, '_', INT, 'methObj75pc.RDS')))
  
  ## Create initial GRange object ##
  ## Reassign column names, avoids clashes with genomation later on ##
  unitecov$original_pos <-  unitecov$start
  unitecov$original_chrom <- unitecov$chr
  colnames(unitecov)
  
  'SLK063_ragtag_chr0.12064403'
  'SLK063_ragtag_chr0 7614 7614 SLK063_ragtag_chr0'
  
  ## Vector containing chromosome, start and end location of CpG ##
  myPos = paste(
    str_extract(E$pos, '^[^\\.]+'),  # Everything before the first full stop
    str_extract(E$pos, '(?<=\\.)[^\\.]+$'),  # Everything after the last full stop
    str_extract(E$pos, '(?<=\\.)[^\\.]+$'),  # Repeated extraction
    str_extract(E$pos, '^[^\\.]+')  # Repeated extraction
  )
  ## Change this vector into a dataframe ##
  df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                   start=sapply(strsplit(myPos, " "), `[`, 2),
                   end=sapply(strsplit(myPos, " "), `[`, 2),
                   original_pos=sapply(strsplit(myPos, " "), `[`, 3),
                   original_chrom=sapply(strsplit(myPos, " "), `[`, 4))
  
  ## Convert this to a GRangesOBJ of every CpG position in the original object ##
  GRangeOBJ <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  
  ## Subset uniteCov object using DMS GRanges object ##
  DMSuniteCov <- selectByOverlap(unitecov, GRangeOBJ)
  
  ## Create percent methylation object ##
  DMSpercMeth <- as.data.frame(percMethylation(DMSuniteCov, rowids = TRUE))
  
  print('##### Save objects #####')
  ## DMS stats dataframe ##
  saveRDS(E, file.path(paste0(OUTDIR, '/', TP, '_', INT, '_PQLseq_DMS_stats.RDS')))
  
  ## DMS percent methylation dataframe ##
  saveRDS(DMSpercMeth, file.path(paste0(OUTDIR, '/', TP, '_', INT, '_PQLseq_DMS_percMeth.RDS')))
  
  ## DMS GRanges object ##
  saveRDS(GRangeOBJ, file.path(paste0(OUTDIR, '/', TP, '_', INT, '_PQLseq_DMS_grange.RDS')))
  
  return(list('stats' = E, 'percMeth' = DMSpercMeth, 'GRangeOBJ' = GRangeOBJ))
}


testRun <- preparePQLseqOutput(TP = 'TP1', INT = 'DS')


TP1DS_summary <- testRun$stats

#### ####
#### TEST CENTER ####
## Quick little PCA, just a little chill PCA guy ##
# Load required package
library(ggplot2)

# Load metadata #
metadata <- data.frame(
  treatment = DMSuniteCov@treatment,
  sampleID  = DMSuniteCov@sample.ids
) %>%
  mutate(
    depth = ifelse(
      treatment == 4, 'shallow', 'deep'
    )
  )

# merge metadata with DMS #
DMSpercMeth2 <- as.data.frame(t(DMSpercMeth)) %>%
  rownames_to_column('sampleID') %>%
  left_join(., metadata, by = 'sampleID') %>%
  column_to_rownames('sampleID') %>%
  dplyr::select(-c(treatment))

# Impute Na values based on TREATMENT MEAN #
df_imputed <- DMSpercMeth2 %>%
  group_by(depth) %>%  # Replace 'Treatment' with your actual column name
  mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  ungroup()

rownames(df_imputed) <- rownames(DMSpercMeth2)



# Perform PCA
pca_result <- prcomp(df_imputed[,1:234], center = TRUE, scale. = TRUE)

# Get variance explained
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create a PCA plot
pca_df <- as.data.frame(pca_result$x)
pca_df$sampleID <- rownames(df_imputed)
pca_df <- pca_df %>%
  left_join(., metadata, by = 'sampleID')

ggplot(pca_df, aes(x = PC1, y = PC2, label = sampleID, color = depth)) +
  geom_point() +
  geom_text(vjust = -1, size = 3) +
  labs(title = "PCA of CpG % Methylation",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_classic()


plot(type = 'bar',
     variance_explained)


#$ -l highmem

