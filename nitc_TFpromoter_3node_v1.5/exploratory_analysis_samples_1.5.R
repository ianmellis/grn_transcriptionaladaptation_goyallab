library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)

datadir <- '~/code/grn_nitc/nitc_TFpromoter_3node_v1.5/'

setwd(datadir)

paramsets <- 1:5

for (paramset in paramsets){

  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T)) %>%
    mutate(paramset = paramset)
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))%>%
    mutate(paramset = paramset)
  
  spstats <- species %>%
    filter((time + 1) %% 300 == 0, # remove when loading subsample
           (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance') %>%
    group_by(mutated_alleles, product) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T))
}
