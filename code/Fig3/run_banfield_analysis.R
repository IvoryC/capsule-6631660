

library(tidyverse)
library(SCRuB)

rescale <- function(m) m / sum(m)
get_rescaled_mat <- function(Q) ( Q %>% apply(MARGIN=1, rescale) %>% t() )


abunds <- read.csv('../data/Banfield/contam_for_liat/YueClareLou_NIHData4Contam_GenomeRelativeAbundanceTable.csv')

ab_mat <- abunds %>% select(sample, genome_name, reads_unfiltered_pairs) %>%
  group_by(sample, genome_name) %>%
  summarise(reads_unfiltered_pairs=sum(reads_unfiltered_pairs)) %>%
  ungroup() %>%
  pivot_wider(id_cols=sample, names_from=genome_name, values_from=reads_unfiltered_pairs)

ab_mat[is.na(ab_mat)] <- 0

ab_mat %<>% remove_rownames %>% column_to_rownames(var="sample")

long_abunds <- ab_mat %>% get_rescaled_mat %>% data.frame() %>%
  mutate(smp_name=row.names(ab_mat)) %>%
  pivot_longer(colnames(ab_mat)[1]:(colnames(ab_mat)[ncol(ab_mat)])) %>%
  filter(value>0)





control_leaks <- c()

for(plate in c(2,3, 4, 5)){#}, 4)){#},4,5)){#} c(2:5)){
  
  p3_locs <- read.csv(paste0('../data/Banfield/contam_for_liat/YueClareLou_NIHData4Contam_PlateLocation/YueClareLou_NIHData4Contam_P', plate, 'SampleMap.csv'))
  p3_locs$cnames <- p3_locs %>% row.names()
  p3_locs %<>% pivot_longer(cols = A:colnames(p3_locs)[ncol(p3_locs)-1], values_to='sample') %>%
    mutate(plate_loc= paste0(name, cnames)) %>%
    select(sample, plate_loc)
  
  metadata <- abunds %>%
    select(sample, sample_type) %>%
    unique() %>%
    mutate(is_control = sample_type=='control') %>%
    merge(p3_locs)
  row.names(metadata) <- metadata$sample
  
  metadata %<>% select(is_control, sample_type, plate_loc)
  
  metadata <- metadata[row.names(metadata)[row.names(metadata) %in% row.names(ab_mat)], ]
  
  
  metadata %>% filter(plate_loc %>% str_detect('L') ) %>% row.names() %>% print
  
  scr_out <- SCRuB( ab_mat[row.names(metadata), ] %>% as.matrix, 
                    metadata %>% mutate(plate_loc=plate_loc %>% str_replace('L', 'Z')), 
                    dist_threshold = 9
  )
  
  print(scr_out$inner_iterations$control$alpha)
  
  
  control_leaks <- c(control_leaks, scr_out$inner_iterations$control$alpha[, scr_out$inner_iterations$control$alpha %>% ncol])
}


strain_comps <- read.csv('../data/Banfield/contam_for_liat/YueClareLou_NIHData4Contam_SamplePairStrainComparisonTable_v2_IncludeNegControls.csv')  

abund_strain_comps <- strain_comps %>% 
  filter(infant1!=infant2) %>%
  merge(long_abunds,
        by.x=c('sample2', 'genome_name'), 
        by.y=c('smp_name', 'name')) %>%
  rbind(  strain_comps %>%  filter(infant1!=infant2) %>% mutate(sample1=sample2,
                                                                sample2=sample1) %>% merge(long_abunds,
                                                                                           by.x=c('sample2', 'genome_name'),
                                                                                           by.y=c('smp_name', 'name'))  
  ) %>% 
  filter(sample1 != sample2) %>%
  filter(popANI>=.99995) %>%
  group_by(sample1, sample2) %>% 
  summarize(n_shared_strains = n(), 
            perc_abund_strain_overlap=sum(value) ) %>%
  ungroup()



nonzero_clonal_controls <- abund_strain_comps %>% filter(sample2 %>% str_detect('L3_EC_001P1'), 
                              sample1 %in% c("L3_082_243G1",
                                             "L3_082_361G1",
                                             "L3_083_038G1" ,
                                             "L3_083_059G1",
                                             "L3_083_088G1" ,
                                             "L3_083_119G1",
                                             "L3_083_239G1",
                                             "L3_EC_001P1") ) %>%
  rbind(  abund_strain_comps %>% filter(sample2 %>% str_detect('L3_EC_002P2'), 
                                        sample1 %in% c('L3_128_124G1', 'L3_128_015G1', 'L3_128_029G1', 'L3_128_020G1', 
                                                       "L3_128_020G1",
                                                       "L3_128_029G1" ,
                                                       "L3_128_037G1",
                                                       "L3_128_043G1",
                                                       "L3_128_054G1",
                                                       "L3_128_070G1",
                                                       "L3_128_124G1",
                                                       "L3_EC_002P2" ) ) ) %>% 
    group_by(sample2) %>% summarize(perc_abund_strain_overlap = max(perc_abund_strain_overlap)) 



data.frame( 'Well-to-well leakage inferred by SCRuB' = 1 - control_leaks, 
            'Fraction of control clonal to most similar sample'= c( 0, nonzero_clonal_controls %>% pull(perc_abund_strain_overlap), 0 ),
            row.names = paste('Plate', c(2:5)) ) %>% t() %>%
  write.csv('../results/Supplementary_figures/F_S8/F_S8_C.csv')








