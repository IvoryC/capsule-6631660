
#### notes ####
# This draft started by taking the first 334 lines of Run_plastma_deconts_and_preds.R

# I thought I would need to run SCRuB for each contamination source, 
# but it looks like SCRuB internalizes that iteration process.

# Note: This ignores plate and plate location to instead handle all samples in one go.

#### libraries ####

# install.packages('splitstackshape', repos='http://cran.us.r-project.org')
# library(splitstackshape)
library(tidyverse)
# library(biomformat) 
library(vegan)
library(glmnet)
library(torch)
install_torch()
library(SCRuB)

sessionInfo()

#### read args ####

# defaults
numTrials = 2 # the number of times to run scrub and produce a result
dropBlanks = .5 # 0-1, the proportion of blanks to drop

# command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
    args2 = args %>% strsplit(split="=")
    key = sapply(args2, function(x) x[[1]])
    argVals = sapply(args2, function(x) x[[2]])
    names(argVals) = key
    if ("trials" %in% key) numTrials = argVals["trials"]
    if ("dropBlanks" %in% key) dropBlanks = argVals["dropBlanks"]
}
message("Will run ", numTrials, " trials.")
message("In each trial, will subsample to remove proportion of blanks: ", dropBlanks)

#### read data ####

full_df = read.csv("../results/data/ivory/before-decontamination/full_df.csv", row.names = 1)
print("full_df.csv has dimentions:")
dim(full_df)

metadata <- read.csv('../data/Fig4_plasma/47212_47212_analysis_mapping.txt', sep='\t')
row.names(metadata) = metadata$X.SampleID
print("metadata has dimentions:")
dim(metadata)

# These lines were not needed when using R version 4.2.2 (2022-10-31), 
# but without them I hit errors when using R version 3.6.3 (2020-02-29)
metadata$sample_plate = as.character(metadata$sample_plate)
metadata$sample_type = as.character(metadata$sample_type)
metadata$sample_well = as.character(metadata$sample_well)

sample_intersect = intersect(row.names(metadata), 
                             row.names(full_df))

counts = full_df[sample_intersect, ]
meta = metadata[sample_intersect, ]

mainResultsDir = "../results/data/ivory/decontamination"
message("Results will be saved under: ", mainResultsDir)

#### run through SCRuB ####

for (trialNumber in 1:numTrials){
    message("Trial ", trialNumber, " of ", numTrials)
    
    seed = trialNumber * 29
    message(paste('Using seed:', seed))
    set.seed(seed)
    
    control_types = c('control blank library prep', 'control blank DNA extraction')
    meta$is_control = meta$sample_type %in% control_types
    
    
    if (dropBlanks > 0){
        totalBlanks = sum(meta$is_control)
        # use ceiling, so if there is at least one blank in the data, 
        # you are guaranteed to have at least one blank to process with.
        numKeepBlanks = ceiling(totalBlanks * (1 - dropBlanks))
        message("Onlye use ", numKeepBlanks, " blanks (out of ", totalBlanks, ")." )
        # split out the controls and sub-sample them, then put back with samples 
        metaSp = split(meta, f=meta$is_control)
        metaSp[['TRUE']] = metaSp[['TRUE']][ sample(1:totalBlanks, size=numKeepBlanks, replace=F), ]
        meta = rbind(metaSp[[1]], metaSp[[2]])
        # subset counts to match meta
        counts = counts[row.names(meta), ]
    }else{
        message("Keeping all blanks.")
    }
    
    
    # when I passed a table with no controls, got error: "Error in rowSums(.) : 'x' must be numeric"
    numControls = sum(meta$is_control)
    message("Data has controls: ", numControls)
    
    scrub_output = SCRuB::SCRuB(counts, metadata = meta[,c("is_control", "sample_type")])
    
    # decontaminated data
    scrub_df <- scrub_output$decontaminated_samples
    
    message("scrub_df has dimensions:")
    dim(scrub_df)
    
    message("Filtering down to columns with sums > 500")
    scrub_df <- scrub_df[, (colSums(scrub_df) > 500) %>% which]
    
    message("scrub_df has dimensions:")
    dim(scrub_df)
    
    numBlanks = sum(meta$is_control)
    resDir = file.path(mainResultsDir, paste0("SCRuB_noPlate_", numBlanks, "-blanks"), paste0("trial_", trialNumber))
    #resDir = paste0("../results/data/ivory/decontamination/SCRuB_noPlate_", numBlanks, "-blanks/trial_", trialNumber)
    suppressWarnings({dir.create(resDir, recursive = T)})
    file = file.path(resDir, "scrub_decontaminated.csv")
    message("Saving scrubbed data as: ", file)
    write.csv(scrub_df, file)
    
}

message("Done!")
