
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
    dir.create("../results/data/ivory/decontamination/raw")
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

#### decontamination step ####

# nope

#### general process ####
# just to match any filtering and/or formating that was done in the decontamination scripts.

message("counts has dimensions:")
dim(counts)

message("Filtering down to columns with sums > 500")
counts_df <- counts[, (colSums(counts) > 500) %>% which]

message("scrub_df has dimensions:")
dim(counts_df)

file = "../results/data/ivory/decontamination/raw/not_decontaminated.csv"
message("Saving scrubbed data as: ", file)
write.csv(counts_df, file)


message("Done!")
