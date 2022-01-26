library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)

datadir <- '~/code/grn_nitc/nitc_TFpromoter_3node_v1.5/'
plotdir <- paste0(datadir, 'exploratory_analysis/')
setwd(datadir)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}

paramsets <- 1:5
allstats <- list()
for (paramset in paramsets){

  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T)) %>%
    mutate(paramset = paramset)
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'.csv'), header = T))%>% # add _q300 when loading subsample
    mutate(paramset = paramset,
           time = 1:nrow(species))
  
  species_sample <- species %>%
    filter((time + 1) %% 300 == 0, # remove when loading subsample
           (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() + 
    facet_grid(mutated_alleles~product)
  
  ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_paramset_', as.character(paramset), '.pdf'))
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T))
  
  if(is.null(dim(allstats))) {
    allstats <- spstats
  } else {
    allstats %<>% bind_rows(spstats)
  }
  
}

write.csv(allstats, file = paste0(plotdir, 'summarystats.csv'))
