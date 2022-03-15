library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(diptest)
library(e1071)
library(ggrepel)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')


# edit as needed
datadir2 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2/samples/'
plotdir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2and5/exploratory_analysis/'
setwd(datadir2)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}
paramsets2 <- 1:9900
lhs_sets2 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets2 %<>%
  mutate(paramset = 1:nrow(lhs_sets2))

allstats2 <- list()
allparams2 <- list()
for (paramset in paramsets2){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_v1.6.2_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              HDTpval = dip.test(abundance)$p.value, .groups = 'keep')
  
  if(is.null(dim(allstats2))) {
    allstats2 <- spstats
  } else {
    allstats2 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams2))) {
    allparams2 <- params
  } else {
    allparams2 %<>% bind_rows(params)
  }
  
}

# edit as needed
datadir5 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.5/samples/'
setwd(datadir5)

paramsets5 <- 1:10000
lhs_sets5 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets5 %<>%
  mutate(paramset = 1:nrow(lhs_sets5))

allstats5 <- list()
allparams5 <- list()
for (paramset in paramsets5){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset) %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_v1.6.5_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              HDTpval = dip.test(abundance)$p.value, .groups = 'keep')
  
  if(is.null(dim(allstats5))) {
    allstats5 <- spstats
  } else {
    allstats5 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams5))) {
    allparams5 <- params
  } else {
    allparams5 %<>% bind_rows(params)
  }
  
}

allstats_full <- bind_rows(allstats2 %>% mutate(version = '1.6.2'), allstats5 %>% mutate(version = '1.6.5'))
allparams_full <- bind_rows(allparams2 %>% mutate(paramset = 1:9900, version = '1.6.2'),
                            allparams5 %>% mutate(paramset = 1:10000, version = '1.6.5'))
lhs_sets_full <- bind_rows(lhs_sets2 %>% mutate(version = '1.6.2'), lhs_sets5 %>% mutate(version = '1.6.5'))

pseud = 0.01

compared_stats <- allstats_full %>% 
  group_by(version, paramset, product, mutated_alleles) %>% 
  pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:HDTpval) %>% 
  pivot_wider(names_from = mutated_alleles, values_from = value) %>% 
  mutate(lfc10 = log2((`1`+pseud)/(`0` + pseud)), 
         delta10 = `1`-`0`,
         lfc21 = log2((`2`+pseud)/(`1` + pseud)), 
         delta21 = `2`-`1`,
         lfc20 = log2((`2`+pseud)/(`0` + pseud)), 
         delta20 = `2`-`0`) %>%
  dplyr::select(-c(`0`:`2`)) %>% 
  pivot_longer(names_to = 'compare', values_to = 'diff', cols = lfc10:delta20) %>%
  inner_join(allstats_full %>% 
               dplyr::select(mutated_alleles, product, paramset, version, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset', 'version')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))

# filter to Hill n < 5