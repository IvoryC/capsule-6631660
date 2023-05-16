
# This draft started by taking a chunk from Run_plastma_deconts_and_preds.R

# install.packages('splitstackshape', repos='http://cran.us.r-project.org')
library(splitstackshape)
library(tidyverse)
library(biomformat) 
library(vegan)
library(glmnet)
library(torch)
# library(microDecon)
# install_torch()

## Load packages ##
require(limma)
require(edgeR)
require(dplyr)
require(snm)
require(doMC)
require(tibble)
require(gbm)

require(parallel)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# BiocManager::install("snm")
# BiocManager::install("edgeR")
# BiocManager::install("limma")
require(edgeR)
require(limma)

sessionInfo()

#### read meta data ####

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

sampleSet = row.names(metadataPSMatchedDPQCFiltered)


#### find data files ####

# SCRuB files in ../results/data/decontaminate/SCRuB/trial_*/scrub_output.csv
infiles = dir("../results/data/ivory/decontamination/SCRuB", pattern="csv", full.names = T, recursive = T)

#### define vsnm function ####

vsnm <- function(qcData){

    numCores <- parallel::detectCores()
    registerDoMC(cores=numCores)
    
    qcMetadata <- metadataPSMatchedDPQCFiltered # ADAPT THIS AS NEEDED
    # qcData <- qqq # ADAPT THIS AS NEEDED
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type_consol + #randomized_disease_type_consol + #disease_type_consol +
                                      host_age + # host_age should be numeric
                                      sex, # sex should be a factor
                                  data = qcMetadata)
    
    # Check row dimensions
    dim(covDesignNorm)[1] == dim(qcData)[1]
    
    # The following corrects for column names that are incompatible with downstream processing
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Quantile normalize and plug into voom
    dge <- edgeR::DGEList(counts = counts)
    vdge <<- limma::voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, 
                  normalize.method="quantile")
    
    # List biological and normalization variables in model matrices
    bio.var <- model.matrix(~disease_type_consol, #randomized_disease_type_consol, #disease_type_consol,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~host_age +
                                sex,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    print(dim(adj.var))
    print(dim(bio.var))
    print(dim(t(vdge$E)))
    print(dim(covDesignNorm))
    
    snmDataObjOnly <- snm::snm(raw.dat = vdge$E, 
                          bio.var = bio.var, 
                          adj.var = adj.var, 
                          rm.adj=TRUE,
                          verbose = TRUE,
                          diagnose = TRUE)
    snmData <<- t(snmDataObjOnly$norm.dat)
    
}

#### use it ####


for (infile in infiles){
    message("Processing file: ", infile)
    outfile = sub(".csv", "_vsnm.csv", infile)
    indata = read.csv(infile, row.names = 1)
    normed <- vsnm(indata[sampleSet, ])
    message("Saving file: ", outfile)
    write.csv(normed, outfile)
}



# scrubbed_normalized <- vsnm(scrub_df[row.names(metadataPSMatchedDPQCFiltered), ])
# 
# microdecon_normalized <- vsnm(microdec_df[row.names(metadataPSMatchedDPQCFiltered), ])
# 
# # scrubbed_normalized <- vsnm(scrub_df[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ), ] )
# # 
# raw_inp <- unique_samps[row.names(metadataPSMatchedDPQCFiltered),]
# raw_inp <- raw_inp[, colSums(raw_inp)>500]
# # raw_inp <- group_to_genus(raw_inp)
# raw_normalized <- vsnm(raw_inp)
# 
# dec_normalized <- vsnm(decontammed_data[row.names(metadataPSMatchedDPQCFiltered), ])#
# # dec_normalized <- vsnm(decontammed_data[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# 
# dec_standard_normalized <- vsnm(decontammed_data_standard[row.names(metadataPSMatchedDPQCFiltered), ])#
# 
# # dec_standard_normalized <- vsnm(decontammed_data_standard[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# dec_lb_normalized <- vsnm(decontammed_data_low_bm[row.names(metadataPSMatchedDPQCFiltered), ])#
# # dec_lb_normalized <- vsnm(decontammed_data_low_bm[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# 
# 
# restrictive_normalized <-  vsnm( restrictive[row.names(metadataPSMatchedDPQCFiltered) %>% as.character(), ])#

message("Done!")
