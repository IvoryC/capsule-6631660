
# This draft started by taking the first 334 lines of Run_plastma_deconts_and_preds.R

# install.packages('splitstackshape', repos='http://cran.us.r-project.org')
# library(splitstackshape)
library(tidyverse)
# library(biomformat) 
library(vegan)
library(glmnet)
library(torch)
# install_torch()

seed=3
## Do several independent iterations. For now, just do it all once.
# for(seed in 0:2){

set.seed(seed)
message(paste('running using seed:', seed))

suppressWarnings({
    dir.create("../results")
    dir.create("../results/data")
    dir.create("../results/data/ivory")
    dir.create("../results/data/ivory/decontamination")
    dir.create("../results/data/ivory/decontamination/SCRuB")
})

#### read data ####

unique_samps <- read.csv("../results/data/ivory/before-decontamination/unique_samps.csv", row.names = 1)
message("unique_samps.csv has dimentions:")
dim(unique_samps)

# metadata must include columns: "sample_type", "sample_well"
unique_metadata <- read.csv("../results/data/ivory/before-decontamination/unique_metadata.csv", row.names = 1)
message("unique_metadata.csv has dimentions:")
dim(unique_metadata)

# Actually.... we just get this from the meta data.
# well_dists <- read.csv("../results/data/ivory/before-decontamination/well_dists.csv", row.names = 1, check.names = F)
# message("well_dists.csv has dimentions:")
# dim(well_dists)

#### run through SCRuB ####

useRPackage = T
if (useRPackage){
    library(SCRuB)
    message(packageVersion("SCRuB"))
    
    control_types = c('control blank library prep', 'control blank DNA extraction')
    unique_metadata$is_control = unique_metadata$sample_type %in% control_types
    
    plates = unique(unique_metadata$sample_plate)
    message("SCRuB the ", length(plates), " plates separately.")
    
    scrub_out_list = list()
    
    for (plate in plates){
        message("SCRuB-ing plate ", plate)
        plate_meta = unique_metadata[unique_metadata$sample_plate == plate,]
        plate_data = unique_samps[row.names(unique_samps) %in% row.names(plate_meta),]
        
        # limit metadata to the EXACT three columns permitted.
        plate_meta = plate_meta[,c("is_control", "sample_type", "sample_well")]
        
        
        # when I passed a table with no controls, got error: "Error in rowSums(.) : 'x' must be numeric"
        # so skip those plates
        numControls = sum(plate_meta$is_control)
        message("Plate has controls: ", numControls)
        message(paste(row.names(plate_meta)[plate_meta$is_control], collapse=", "))
        
        if (numControls > 0 ){
            scrub_output = SCRuB::SCRuB(plate_data, metadata = plate_meta)
            
            # decontaminated data
            scrub_out_list[[plate]] <- scrub_output$decontaminated_samples

            # I don't understand this error: "Error: value for ‘self’ not found"
            # happened to plate CMI_Plasma_12692_P4 -- but not the following time. (?)
            # only happens some of the time.
        }else{
            message("Skipping plate")
            scrub_out_list[[plate]] <- plate_data
        }
    }
    scrub_df = do.call("rbind", scrub_out_list)
    
}else{
    source('SCRuB/lsq_initializations.R')
    source('SCRuB/spatial_functions.R')
    source('SCRuB/main_functions.R')
    
    well_dists <- unique_metadata %>%
        mutate(well_loc= sample_well %>% substr(1,1) %>% sapply( function(x) which(LETTERS==x)[1]), 
               indices_loc = sample_well%>% substr(2,3) %>% as.integer ) %>%
        select(well_loc, indices_loc) %>%
        dist(method = 'euclidean') %>% as.matrix()
    
    well_dists = round(well_dists, digits = 2)
    
    
    tmp_fwd <- as.matrix( unique_samps[ unique_metadata %>% 
                                            filter(sample_type %>% 
                                                       str_detect('control blank library prep') == F) %>% 
                                            pull(X.SampleID) %>% 
                                            as.character(), ] )
    
    ######## encounters problem here
    plasma_scrubbed <- spatial_SCRUB(data=tmp_fwd, 
                                     is_control = unique_metadata[row.names(tmp_fwd),]$sample_type=='control blank DNA extraction',
                                     well_dists = well_dists, 
                                     dist_threshold =1.5
    )
    
    Qrow.names( plasma_scrubbed$decontaminated_samples ) <- unique_metadata[row.names(tmp_fwd),
    ][ unique_metadata[row.names(tmp_fwd),
    ]$sample_type!='control blank DNA extraction', ] %>% row.names()
    
    
    
    next_lvl_mat <- plasma_scrubbed$decontaminated_samples %>%
        rbind(unique_samps[row.names(unique_metadata %>%
                                         filter(sample_type=='control blank library prep')), ])
    
    
    
    plasma_scrubbed_2nd <- spatial_SCRUB(data=next_lvl_mat, 
                                         is_control = c( rep(F, nrow(plasma_scrubbed$decontaminated_samples) ),
                                                         rep(T, sum(unique_metadata$sample_type=='control blank library prep'))),
                                         well_dists = well_dists,
                                         dist_threshold =  1.5 )
    
    
    row.names(plasma_scrubbed_2nd$decontaminated_samples) <- unique_metadata[row.names(plasma_scrubbed$decontaminated_samples),
    ][ unique_metadata[row.names(plasma_scrubbed$decontaminated_samples ),
    ]$sample_type!='control blank library prep', ] %>% row.names()
    
    
    colnames(plasma_scrubbed_2nd$decontaminated_samples) <- colnames(unique_samps)
    
    scrub_df <- plasma_scrubbed_2nd$decontaminated_samples[,  ( ( plasma_scrubbed_2nd$decontaminated_samples %>% colSums() ) > 0 ) %>% which]
    scrub_df <- scrub_df[row.names(metadataPSMatchedDPQCFiltered), ]
    
    
}



message("scrub_df has dimensions:")
dim(scrub_df)

message("Filtering down to columns with sums > 500")
scrub_df <- scrub_df[, (colSums(scrub_df) > 500) %>% which]

message("scrub_df has dimensions:")
dim(scrub_df)

# trialNumber=seed
# file=paste0("../results/data/ivory/decontamination/SCRuB/trial_", trialNumber, "/scrub_output.csv")
file="../results/data/ivory/decontamination/SCRuB/scrub_output.csv"
message("Saving scrubbed data as: ", file)
write.csv(scrub_df, file)
