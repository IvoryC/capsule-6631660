
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
    dir.create("../results/data/ivory/decontamination/SCRuB")
})

#### read data ####

# unique_samps <- read.csv("../results/data/ivory/before-decontamination/unique_samps.csv", row.names = 1)
# message("unique_samps.csv has dimentions:")
# dim(unique_samps)

full_df = read.csv("../results/data/ivory/before-decontamination/full_df.csv", row.names = 1)
message("full_df.csv has dimentions:")
dim(full_df)

## metadata must include columns: "sample_type", "sample_well", "sample_plate"
# unique_metadata <- read.csv("../results/data/ivory/before-decontamination/unique_metadata.csv", row.names = 1)
# message("unique_metadata.csv has dimentions:")
# dim(unique_metadata)
metadata <- read.csv('../data/Fig4_plasma/47212_47212_analysis_mapping.txt', sep='\t')
row.names(metadata) = metadata$X.SampleID
message("metadata has dimentions:")
dim(metadata)

sample_intersect = intersect(row.names(metadata), 
                             row.names(full_df))

full_df = full_df[sample_intersect, ]
metadata = metadata[sample_intersect, ]

# These lines were not needed when using R version 4.2.2 (2022-10-31), 
# but without them I hit errors when using R version 3.6.3 (2020-02-29)
# unique_metadata$sample_plate = as.character(unique_metadata$sample_plate)
# unique_metadata$sample_type = as.character(unique_metadata$sample_type)
# unique_metadata$sample_well = as.character(unique_metadata$sample_well)

# Actually.... we just get this from the meta data.
# well_dists <- read.csv("../results/data/ivory/before-decontamination/well_dists.csv", row.names = 1, check.names = F)
# message("well_dists.csv has dimentions:")
# dim(well_dists)

#### run through SCRuB ####

numTrials = 10
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
    metadata$is_control = metadata$sample_type %in% control_types
    
    # Note: This decontaminates per plate. And one plate, is not decontaminated all. But the well location info is optional.  
    #       Could I instead decontaminate the whole data set as one if I just ignore the well location info?
    
    plate_meta_split = split(metadata, f=metadata$sample_plate)
    message("SCRuB the ", length(plate_meta_split), " plates separately.")
    
    scrub_out_list = list()
    
    for (plate in names(plate_meta_split)){
        message("SCRuB-ing plate ", plate)
        plate_meta = plate_meta_split[[plate]]
        plate_data = full_df[row.names(plate_meta),]
        
        # limit metadata to the EXACT three columns permitted.
        plate_meta = plate_meta[,c("is_control", "sample_type", "sample_well")]
        
        
        # when I passed a table with no controls, got error: "Error in rowSums(.) : 'x' must be numeric"
        # so skip any plate with 0 controls
        numControls = sum(plate_meta$is_control)
        message("Plate has controls: ", numControls)
        message(paste(row.names(plate_meta)[plate_meta$is_control], collapse=", "))
        
        if (numControls > 0 ){
            scrub_output = SCRuB::SCRuB(plate_data, metadata = plate_meta)
            
            # decontaminated data
            scrub_out_list[[plate]] <- scrub_output$decontaminated_samples
            
        }else{
            message("Skipping plate")
            scrub_out_list[[plate]] <- plate_data
        }
    }
    scrub_df = do.call("rbind", scrub_out_list)
    
    message("scrub_df has dimensions:")
    dim(scrub_df)
    
    message("Filtering down to columns with sums > 500")
    scrub_df <- scrub_df[, (colSums(scrub_df) > 500) %>% which]
    
    message("scrub_df has dimensions:")
    dim(scrub_df)
    
    resDir = paste0("../results/data/ivory/decontamination/SCRuB/trial_", trialNumber)
    suppressWarnings({dir.create(resDir)})
    file = file.path(resDir, "scrub_output.csv")
    message("Saving scrubbed data as: ", file)
    write.csv(scrub_df, file)
    
}

message("Done!")
