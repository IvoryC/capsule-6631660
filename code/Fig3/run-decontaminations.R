library(SCRuB)
rescale <- function(x) x/sum(x)
torch::install_torch()

library(decontam)
library(microDecon)
CLEAN_SAMPLES_DECONTAM<- function(samples, 
                                  controls, 
                                  threshold = 0.1){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  
  ps <- samples %>% rbind(controls)
  
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isContaminant(ps, neg=neg, threshold = threshold,
                       detailed = T,
                       method='prevalence'
  )
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == T] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$contaminant),
              cont_idx = out, 
              estimated_sources = controls, 
              gamma=1 - out$p ) )
}

CLEAN_SAMPLES_DECONTAM_low_biomass <- function(samples, 
                                               controls, 
                                               threshold = 0.5){
  n_sinks <- nrow(samples)
  n_contams <- nrow(controls)
  
  ps <- samples %>% rbind(controls)
  
  neg <-c(rep(FALSE, n_sinks), rep(TRUE, n_contams))
  
  out <- isNotContaminant(ps, neg, threshold = threshold, 
                          detailed = T)
  
  remove_out_samps <- function(mat, out){
    mat %>% apply(MARGIN = 1, function(x){
      x[out == FALSE] <- 0 
      return(x)
    }  ) %>% t()
  }
  return(list(decontaminated_samples = remove_out_samps(samples, out$not.contaminant),
              cont_idx = out, 
              estimated_sources = controls, 
              gamma = out$p) )
}

CLEAN_SAMPLES_MICRODECON <- function(smps, cnts){
  tmp <- data.frame( rbind( colnames(smps), cnts, smps ) %>% t() )
  tmp[, 2:( ncol(tmp) ) ]  <- as.numeric( as.matrix( tmp[, 2:( ncol(tmp) ) ] ) )
  decontaminated <- decon(data = tmp, numb.blanks=nrow(cnts), numb.ind= c(nrow(smps)), taxa=F)
  
  md_out <- smps*0
  md_out[,decontaminated$decon.table$V1 %>% as.character() %>% unname()] <- decontaminated$decon.table[,3:(2 + nrow(smps) )] %>% t()
  
  return(list(decontaminated_samples=md_out, 
              gamma = (1 - md_out /  smps ) %>% colMeans(na.rm = T) 
  ))
}


REMOVE_ALL_CONTROL_SPECS <- function(samples, controls){
  X = samples
  to_remove <- which( (controls %>% colSums() )  > 0  )
  X[,to_remove] <- 0
  return(list(decontaminated_samples=X, 
              gamma=colSums(controls) > 0))
}

alternative_approach_funcs <- list('Decontam'=CLEAN_SAMPLES_DECONTAM, 
                                   'Decontam (LB)'=CLEAN_SAMPLES_DECONTAM_low_biomass,
                                   'microDecon'=CLEAN_SAMPLES_MICRODECON,
                                   'Restrictive'=REMOVE_ALL_CONTROL_SPECS
)

extraction_metadata1 <- read.csv('../data/Fig3/extraction_metadata.csv', row.names=1) %>% mutate(is_control = as.logical(is_control), 
                                                                                                 well_loc=as.character(well_loc))
extraction_metadata1['neg2_9_S56', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'G12')
extraction_metadata1['neg2_10_S47', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'H12')

extraction_metadata2 <- read.csv('../data/Fig3/extraction_metadata2.csv', row.names=1) %>% mutate(is_control = as.logical(is_control), 
                                                                                                  well_loc=as.character(well_loc))

extraction_metadata2['neg2_3_2_S50', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'F2')

extraction_metadata2['zneg2_9_S100', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'G12')
extraction_metadata2['zneg2_10_S93', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'H12')


dilution_samples <- read.csv('../data/Fig3/dilution-samples.csv', row.names=1) %>% as.matrix
colnames(dilution_samples) <- read.csv('../data/Fig3/20220926_16S_Decontamination_MK.220929.asvTable.csv') %>% mutate(Genus=Genus) %>% pull(Genus)


metadata <- read.csv('../data/Fig3/metadata.csv', row.names=1) %>%
  mutate(is_control=is_control %>% as.logical , 
         well_loc=as.character(well_loc))
data <- read.csv('../data/Fig3/ASV-table.csv', row.names=1) %>% as.matrix

geni_cols <-  read.csv('../data/Fig3/20220926_16S_Decontamination.220929.asvTable.csv') %>% mutate(Genus=Genus) %>% pull(Genus)


data <- data[ unique( c( data[ rowSums(data) > 5000, ] %>% row.names(),  
                         data[ str_detect(row.names(data), 'neg' ), ] %>% row.names() ) ), ]

data <- data[rowSums(data) > 500,  ]

metadata <- metadata[row.names(data), ]



extraction_metadata2 <- extraction_metadata2[row.names(extraction_metadata2) %in% row.names(data), ]
print(dim(extraction_metadata2))

set.seed(2022)
inds <- which( (metadata$Groups == 'Lib_neg')|(metadata$is_control==F) )
scr_lib_out <- SCRuB(data[inds, ] , 
                     metadata[inds, 1:3], 
                     control_order= c('Lib_neg'), 
                     dist_threshold = 1.5 
)


extraction_metadata1 <- extraction_metadata1[row.names(extraction_metadata1) %in% row.names(data), ]
next_level_setup <- ( scr_lib_out$decontaminated_samples %>% rbind(
  data[ row.names(data)[ ( row.names( data ) %in% row.names(scr_lib_out$decontaminated_samples) == F )&
                           ( row.names( data ) %in% row.names( extraction_metadata1  ) ) ], ]
) )[row.names( extraction_metadata1 ), ]

extraction_out_1 <- SCRuB(next_level_setup, 
                          extraction_metadata1[, 1:3])
             

gamma_df <- data.frame( Genus=geni_cols, 
                        gamma=extraction_out_1$inner_iterations$neg$gamma,
                        neg_abunds=data[ data %>% row.names() %>% str_starts('neg'), ] %>% colSums(), 
                        samp_abunds=data[extraction_out_1$decontaminated_samples %>% row.names, ] %>% colSums(), 
                        n_samp_presences= ( data[extraction_out_1$decontaminated_samples %>% row.names, ] > 0 ) %>% colSums()
) %>%
  group_by(Genus) %>% summarize(
    neg_abunds=sum(neg_abunds), 
    samp_abunds=sum(samp_abunds),
    n_samp_presences=max(n_samp_presences),
    SCRuB=sum(gamma)
  ) %>% 
  ungroup() %>%
  merge( 
    data.frame( Genus=colnames(dilution_samples), zymo=dilution_samples['MK1_0_S5', ] ) %>%
      group_by(Genus) %>% 
      summarize(zymo=sum(zymo)) %>% 
      ungroup(), 
    by='Genus'
  )


tmp_data <- data[row.names(extraction_metadata1), ]

for(func in names(alternative_approach_funcs)){
  
  gamma_df <- gamma_df %>% merge(
    data.frame( Genus=geni_cols, 
                func = replace_na( alternative_approach_funcs[[func]](tmp_data[extraction_metadata1$is_control==F, ],
                                                                      tmp_data[extraction_metadata1$is_control==T, ])$gamma, 0) ) %>%
      group_by(Genus) %>% 
      summarize(func=mean(func)) %>% 
      ungroup() %>%
      mutate(!!func:=func) %>% 
      select(-func),
    by='Genus') 
}


gamma_df %>% write_csv('../results/data/Fig_3/predicted_contaminants.csv')



library(pROC)
roc.test( roc(gamma_df %>% mutate(label=zymo>200), label, SCRuB), 
          roc(gamma_df %>% mutate(label=zymo>200), label, Decontam) )

roc.test( roc(gamma_df %>% mutate(label=zymo>200), label, SCRuB), 
          roc(gamma_df %>% mutate(label=zymo>200), label, microDecon) )

roc.test( roc(gamma_df %>% mutate(label=zymo>200), label, SCRuB), 
          roc(gamma_df %>% mutate(label=zymo>200), label, `Decontam (LB)`) )

roc.test( roc(gamma_df %>% mutate(label=zymo>200), label, SCRuB), 
          roc(gamma_df %>% mutate(label=zymo>200), label, Restrictive) )





print(dim(extraction_metadata2))

extraction_metadata2['neg2_3_2_S50', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'F2')

extraction_metadata2['zneg2_9_S100', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'G12')
extraction_metadata2['zneg2_10_S93', c('is_control', 'Groups', 'well_loc')] <-  c(TRUE, 'neg', 'H12')
extraction_metadata2 <- extraction_metadata2[row.names(extraction_metadata2) %in% row.names(data), ]

next_level_setup <- ( scr_lib_out$decontaminated_samples %>% rbind(
  data[ row.names(data)[ ( row.names( data ) %in% row.names(scr_lib_out$decontaminated_samples) == F )&
                           ( row.names( data ) %in% row.names(extraction_metadata2) ) ], ]
) )[row.names(extraction_metadata2), ] 

print(dim(extraction_metadata2))

extraction_out_2 <- SCRuB(next_level_setup, 
                          extraction_metadata2[, 1:3])


rbind(extraction_out_1$decontaminated_samples, extraction_out_2$decontaminated_samples) %>% 
  write.csv('../results/data/Fig_3/decontaminated-samples/SCRuBout.csv')



alternative_approach_funcs <- list('Decontam'=CLEAN_SAMPLES_DECONTAM, 
                                   'Decontam (LB)'=CLEAN_SAMPLES_DECONTAM_low_biomass,
                                   'microDecon'=CLEAN_SAMPLES_MICRODECON,
                                   'Restrictive'=REMOVE_ALL_CONTROL_SPECS
)


inds_2 <- row.names(metadata)

data <- data[inds_2, ]
metadata <- metadata[inds_2, ]


alternate_decontaminations <- alternative_approach_funcs %>%
  lapply( function(func){ func(data[metadata$is_control==F, ],
                               data[metadata$is_control==T, ])$decontaminated_samples })


for( method in c('Decontam', 'Decontam (LB)', 'microDecon', 'Restrictive' )){
  alternate_decontaminations[[method]] %>% write.csv(paste0('../results/data/Fig_3/decontaminated-samples/', method, 'out.csv'))
}

data.frame(extraction_out_2$inner_iterations$neg$gamma) %>% write_csv('../results/data/Fig_3/predicted_plate1_contaminant.csv')
