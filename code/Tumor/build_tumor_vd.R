library(tidyverse)
library(ggvenn)
library(SCRuB)


load('../data/Fig4/cancer_decont_input.Rdata')
# define names of tissues
control_tissues <- c( 'control', 'NTC', 'paraf control' )
sample_tissues <- c( 'Bone', 'Breast','Colon', 'GBM', 'Lung','Melanoma', 'Ovary', 'Pancreas')


their_output <- (read.csv('../data/Fig4/df_freqs.csv') ) %>% t()
their_output[is.na(their_output)] <- 0

tissue_names <- all_info %>% filter(tissue.type %in% sample_tissues) %>% pull(new_SEQ_NAMES)

row.names( their_output ) <- row.names( their_output ) %>% str_replace('X2', '2')
their_output[their_output %>% row.names() %in% pull( all_info, new_SEQ_NAMES), ] %>% nrow()
filter_tracker <- read.csv('../data/Fig4/df_filter_tracker.csv')

their_output <-  their_output[, filter_tracker %>% pull(species)!= '' ]

all_inf_phylos <- read.csv('../data/Fig4/df_freqs_before_global_contaminants_filter.csv')[2:8]

all_inf_phylos$all_inf_idx <- seq(1, nrow(all_inf_phylos))
filter_tracker_phylos <- filter(filter_tracker, species != '')[, 2:8]
filter_tracker_phylos$filt_idx <- seq(1, nrow(filter_tracker_phylos))

merged_idx <- all_inf_phylos %>%
  merge(filter_tracker_phylos, 
        by.x= c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
        by.y= c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
        all=TRUE) %>% 
  select( all_inf_idx, filt_idx ) %>%
  arrange(all_inf_idx)

merged_idx$filt_idx[ merged_idx$filt_idx %>% is.na() ] <- seq( 1 + nrow(filter_tracker_phylos), nrow(all_inf_phylos) )


num_pads <- nrow(all_inf_phylos) - nrow(filter_tracker_phylos)
their_out_padded <- their_output %>% cbind( matrix( 0,
                                                    their_output %>% nrow(), 
                                                    num_pads) )

matching_format <- their_out_padded[2:nrow(their_out_padded), merged_idx$filt_idx]

## build plots
decontam_cont_out_lb <- read.csv( file = '../results/data/Tumor/decontam_low_biomass_result.csv', row.names=1)
reverse_cont_out <- read.csv(file = '../results/data/Tumor/SCRUB_result.csv', row.names=1)
global_fully_restrictive <- read.csv('../results/data/Tumor/global_fully_restrictive_result.csv', row.names=1)


melanoma_idcs <- which( ( ( all_info %>%
                              filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)]  ) %>% as.matrix() %>%
                          colSums() > 1e-2 ) %>% unname()



# dir.create('../results/Plots/')
# dir.create('../results/Plots/fig_3/')

ggvenn(
  list(
    "Custom"= ((matching_format[which( row.names(matching_format) %in% 
                                         ( all_info %>% filter( tissue.type == 'Melanoma') )$new_SEQ_NAMES ), ] %>%
                  colSums )[melanoma_idcs] == 0 ) %>% which , 
    "SCRUB" = ( ( ( ( reverse_cont_out %>%
                        filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                    colSums )[melanoma_idcs] == 0 ) %>% which ,
    "Decontam (LB)"=( ( ( ( decontam_cont_out_lb %>%
                              filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                          colSums )[melanoma_idcs] == 0 ) %>% which, 
    "Restrictive"=( ( ( ( global_fully_restrictive %>%
                            filter(tissue.type=='Melanoma') )[, 6:ncol(reverse_cont_out)] ) %>%
                        colSums )[melanoma_idcs] == 0 ) %>% which
  ), 
  columns = c("Custom", 
              "SCRUB" , 
              "Decontam (LB)", 
              "Restrictive"),
  fill_color=c('#FF69B4', '#c44e52', 'darkgreen', '#dd8452'),
  show_percentage = F, 
  fill_alpha = .7
)
ggsave('../results/Supplementary_figures/F_S10/F_S10_G_with_numbers.pdf', device='pdf', dpi=900)


ggvenn(
  list( "Custom"= 1:10, "SCRUB" =1:10,"Decontam (LB)"=1:10,"Restrictive"=1:10),
  columns = c("Restrictive (Original)" ,
              "SCRUB" ,
              "Decontam (LB)",
              "Restrictive"),
  fill_color=c('#FF69B4', '#c44e52', 'darkgreen', '#dd8452'),
  show_percentage = F, 
  fill_alpha = .8, 
  show_elements=F, 
  text_size=0, 
  set_name_size=0 ) 

ggsave('../results/Supplementary_figures/F_S10/F_S10_G_without_numbers.pdf', device='pdf', dpi=900)







load('../data/Fig4/cancer_decont_input.Rdata')
set.seed(1)
specific_batch_data <- all_info %>% filter(PCR...New.Batch == 2063, #2005, 
                                             tissue.type %in% c('Melanoma', 'NTC'))
row.names(specific_batch_data) <- specific_batch_data$new_SEQ_NAMES
metadata <- specific_batch_data %>% select(new_SEQ_NAMES, tissue.type)
metadata$is_control = metadata$tissue.type=='NTC'
metadata <- metadata %>% select(is_control, tissue.type)

ccr_samples <- specific_batch_data[, 6:ncol(specific_batch_data)] %>% as.matrix
ccr_samples <- ccr_samples[,colSums(ccr_samples)>0]

ccr_samples <- specific_batch_data[, 6:ncol(specific_batch_data)] %>% as.matrix
ccr_samples <- ccr_samples[,colSums(ccr_samples)>0]

scr_out <-  SCRuB(ccr_samples, 
                  metadata = metadata)

most_contaminated <- ccr_samples[ (scr_out$p %>% sort() ) %>% head(4) %>% names(), ]

least_contaminated <-ccr_samples[ (scr_out$p %>% sort() ) %>% tail(4) %>% names(), ]

scr_out$p %>% sort() 
  
ntcs <- ccr_samples[filter(metadata, is_control)%>% row.names, ] 

cont_elements <- which( (ntcs %>% colSums() > 1000)&( (ntcs>0) %>% colSums() >= 2 ) ) %>% names

rbind( most_contaminated, 
       least_contaminated, 
       ntcs) %>%
  write.csv('../results/data/Tumor/selected_melanomas.csv')

