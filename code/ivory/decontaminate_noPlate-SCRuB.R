
# This draft started by taking the first 334 lines of Run_plastma_deconts_and_preds.R

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

suppressWarnings({
    dir.create("../results")
    dir.create("../results/data")
    dir.create("../results/data/ivory")
    dir.create("../results/data/ivory/decontamination")
    dir.create("../results/data/ivory/decontamination/SCRuB_noPlate")
})

#### read data ####

full_df = read.csv("../results/data/ivory/before-decontamination/full_df.csv", row.names = 1)
message("full_df.csv has dimentions:")
dim(full_df)

metadata <- read.csv('../data/Fig4_plasma/47212_47212_analysis_mapping.txt', sep='\t')
row.names(metadata) = metadata$X.SampleID
message("metadata has dimentions:")
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

#### run through SCRuB ####

numTrials = 2
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) numTrials = args[1]

for (trialNumber in 1:numTrials){
    message("Trial ", trialNumber, " of ", numTrials)
    
    seed = trialNumber * 29
    message(paste('Using seed:', seed))
    set.seed(seed)
    
    # I thought I would need to run SCRuB for each contamination source, 
    # but it looks like SCRuB internalizes that iteration process.
    control_types = c('control blank library prep', 'control blank DNA extraction')
    meta$is_control = meta$sample_type %in% control_types
    
    # Note: This ignores plate and plat location to instead handle all samples in one go.
    
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
    
    resDir = paste0("../results/data/ivory/decontamination/SCRuB_noPlate/trial_", trialNumber)
    suppressWarnings({dir.create(resDir)})
    file = file.path(resDir, "scrub_decontaminated.csv")
    message("Saving scrubbed data as: ", file)
    write.csv(scrub_df, file)
    
}

message("Done!")
