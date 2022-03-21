library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(diptest)
library(e1071)
library(ggrepel)
library(corrplot)
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
compared_stats %<>% ungroup() %>% inner_join(lhs_sets_full %>% filter(Hill_coefficient_n < 5) %>% dplyr::select(version, paramset), by = c('version','paramset'))

# plot wt/mut summary stats against mean expression

unistats<-unique(compared_stats$stat)

for (st in unistats) {
  
  pvs1 <- ggplot(allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(mean_product, eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = 10), linetype = 2) +
    # geom_text(aes(mean_product, eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic() +
    ggtitle(paste0(st, ' vs mean, only including Hill n < 5'))
  
  pvs2 <- ggplot(allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n<5)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic()  +
    ggtitle(paste0(st, ' vs mean, only including Hill n < 5'))
  
  pvs3 <- ggplot(allstats_full %>% 
                   inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% 
                   filter(product == 'B1', Hill_coefficient_n<5, mutated_alleles == 1) %>% 
                   mutate(is10 = mean_product > 10)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)), color = is10), alpha = 0.3, stroke = 0) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    scale_color_manual(values = c('grey50', 'black')) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    # facet_grid(~mutated_alleles) +
    theme_classic()  +
    theme(legend.position = 'none') +
    ylab(st) +
    ggtitle(paste0(st, ' vs mean, only including Hill n < 5'))
  
  ggsave(pvs1, file = paste0(plotdir, 'PerGenotype_', st, '_vs_mean.pdf'), width = 16, height = 8) 
  ggsave(pvs2, file = paste0(plotdir, 'PerGenotype_', st, '_vs_logmean.pdf'), width = 16, height = 8) 
  ggsave(pvs3, file = paste0(plotdir, 'WTMUT_', st, '_vs_logmean.pdf'), width = 4, height = 4) 
}

# get example histograms and traces for same-mean classes (high vs low Bimod and high vs low CV)
# aim for mean of roughly 50 (+/- 5)
high_bimod_wtmut <- allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1) %>% arrange(-bimodality_coef) %>% head(10)    
low_bimod_wtmut <- allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1) %>% arrange(bimodality_coef) %>% head(10)

high_cv_wtmut <- allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1) %>% arrange(-cv_product) %>% head(10)  
low_cv_wtmut <- allstats_full %>% inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1) %>% arrange(cv_product) %>% head(10)

high_mean_wtmut <- allstats_full %>% 
  inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% 
  filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 15, mutated_alleles == 1, cv_product > 1.9, cv_product < 2.1) %>% 
  arrange(-mean_product) %>% head(10)  
low_mean_wtmut <- allstats_full %>% 
  inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% 
  filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 15, mutated_alleles == 1, cv_product > 1.9, cv_product < 2.1) %>% 
  arrange(mean_product) %>% head(10)  

if(!dir.exists(paste0(plotdir, 'traces/'))){
  dir.create(paste0(plotdir, 'traces/'))
}

bimod_vs_cv_wtmut_mean50B1 <- ggplot() + 
  geom_point(data = allstats_full %>% 
               inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% 
               filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1), 
             aes(cv_product, bimodality_coef)) +
  geom_text_repel(data = allstats_full %>% 
               inner_join(lhs_sets_full, by = c('version', 'paramset')) %>% 
               filter(product == 'B1', Hill_coefficient_n < 5, mean_product > 45,mean_product < 55, mutated_alleles == 1), 
             aes(cv_product, bimodality_coef, label = paste0(version, '_', paramset))) +
  theme_bw() +
  ggtitle('B1 in WT/MUT genotype, 45 < mean < 55\nCV vs Bimodality coefficient')
ggsave(bimod_vs_cv_wtmut_mean50B1, file = paste0(plotdir, 'bimod_vs_cv_wtmut_mean50B1.pdf'))

for (pind in 1:nrow(high_bimod_wtmut)) {
  
  ver = high_bimod_wtmut$version[pind]
  pset = high_bimod_wtmut$paramset[pind]
  
  if(ver == '1.6.2') {
    tracedir = datadir2
  } else {
    tracedir = datadir5
  }
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  traceplot0b <- plot_traces_ver_Bfocus(species, 7500, 8000, 'WT/WT')
  traceplot1b <- plot_traces_ver_Bfocus(species, 107500, 108000, 'WT/MUT')
  traceplot2b <- plot_traces_ver_Bfocus(species, 207500, 208000, 'MUT/MUT')
  
  f0b<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus.pdf')
  f1b<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus.pdf')
  f2b<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus.pdf')
  
  ggsave(plot = traceplot0b, f0b, width = 8, height = 6)
  ggsave(plot = traceplot1b, f1b, width = 8, height = 6)
  ggsave(plot = traceplot2b, f2b, width = 8, height = 6)
  
  
  traceplot0s <- plot_traces_ver(species, 7500, 7600, 'WT/WT')
  traceplot1s <- plot_traces_ver(species, 107500, 107600, 'WT/MUT')
  traceplot2s <- plot_traces_ver(species, 207500, 207600, 'MUT/MUT')
  
  f0s<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_short.pdf')
  f1s<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_short.pdf')
  f2s<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_short.pdf')
  
  ggsave(plot = traceplot0s, f0s, width = 8, height = 6)
  ggsave(plot = traceplot1s, f1s, width = 8, height = 6)
  ggsave(plot = traceplot2s, f2s, width = 8, height = 6)
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.pdf')
  f1sb<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.pdf')
  f2sb<-paste0(plotdir, 'traces/highBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.pdf')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter(time %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir, 'traces/highBimodality_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir, 'traces/highBimodality_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
}



for (pind in 1:nrow(low_bimod_wtmut)) {
  
  ver = low_bimod_wtmut$version[pind]
  pset = low_bimod_wtmut$paramset[pind]
  
  if(ver == '1.6.2') {
    tracedir = datadir2
  } else {
    tracedir = datadir5
  }
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  traceplot0b <- plot_traces_ver_Bfocus(species, 7500, 8000, 'WT/WT')
  traceplot1b <- plot_traces_ver_Bfocus(species, 107500, 108000, 'WT/MUT')
  traceplot2b <- plot_traces_ver_Bfocus(species, 207500, 208000, 'MUT/MUT')
  
  f0b<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus.pdf')
  f1b<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus.pdf')
  f2b<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus.pdf')
  
  ggsave(plot = traceplot0b, f0b, width = 8, height = 6)
  ggsave(plot = traceplot1b, f1b, width = 8, height = 6)
  ggsave(plot = traceplot2b, f2b, width = 8, height = 6)
  
  
  traceplot0s <- plot_traces_ver(species, 7500, 7600, 'WT/WT')
  traceplot1s <- plot_traces_ver(species, 107500, 107600, 'WT/MUT')
  traceplot2s <- plot_traces_ver(species, 207500, 207600, 'MUT/MUT')
  
  f0s<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_short.pdf')
  f1s<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_short.pdf')
  f2s<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_short.pdf')
  
  ggsave(plot = traceplot0s, f0s, width = 8, height = 6)
  ggsave(plot = traceplot1s, f1s, width = 8, height = 6)
  ggsave(plot = traceplot2s, f2s, width = 8, height = 6)
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.pdf')
  f1sb<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.pdf')
  f2sb<-paste0(plotdir, 'traces/lowBimodality_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.pdf')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter(time %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir, 'traces/lowBimodality_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir, 'traces/lowBimodality_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
}


for (pind in 1:nrow(high_cv_wtmut)) {
  
  ver = high_cv_wtmut$version[pind]
  pset = high_cv_wtmut$paramset[pind]
  
  if(ver == '1.6.2') {
    tracedir = datadir2
  } else {
    tracedir = datadir5
  }
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  traceplot0b <- plot_traces_ver_Bfocus(species, 7500, 8000, 'WT/WT')
  traceplot1b <- plot_traces_ver_Bfocus(species, 107500, 108000, 'WT/MUT')
  traceplot2b <- plot_traces_ver_Bfocus(species, 207500, 208000, 'MUT/MUT')
  
  f0b<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus.pdf')
  f1b<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus.pdf')
  f2b<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus.pdf')
  
  ggsave(plot = traceplot0b, f0b, width = 8, height = 6)
  ggsave(plot = traceplot1b, f1b, width = 8, height = 6)
  ggsave(plot = traceplot2b, f2b, width = 8, height = 6)
  
  
  traceplot0s <- plot_traces_ver(species, 7500, 7600, 'WT/WT')
  traceplot1s <- plot_traces_ver(species, 107500, 107600, 'WT/MUT')
  traceplot2s <- plot_traces_ver(species, 207500, 207600, 'MUT/MUT')
  
  f0s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_short.pdf')
  f1s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_short.pdf')
  f2s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_short.pdf')
  
  ggsave(plot = traceplot0s, f0s, width = 8, height = 6)
  ggsave(plot = traceplot1s, f1s, width = 8, height = 6)
  ggsave(plot = traceplot2s, f2s, width = 8, height = 6)
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.pdf')
  f1sb<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.pdf')
  f2sb<-paste0(plotdir, 'traces/highCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.pdf')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter(time %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir, 'traces/highCV_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir, 'traces/highCV_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
}



for (pind in 1:nrow(low_cv_wtmut)) {
  
  ver = low_cv_wtmut$version[pind]
  pset = low_cv_wtmut$paramset[pind]
  
  if(ver == '1.6.2') {
    tracedir = datadir2
  } else {
    tracedir = datadir5
  }
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  traceplot0b <- plot_traces_ver_Bfocus(species, 7500, 8000, 'WT/WT')
  traceplot1b <- plot_traces_ver_Bfocus(species, 107500, 108000, 'WT/MUT')
  traceplot2b <- plot_traces_ver_Bfocus(species, 207500, 208000, 'MUT/MUT')
  
  f0b<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus.pdf')
  f1b<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus.pdf')
  f2b<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus.pdf')
  
  ggsave(plot = traceplot0b, f0b, width = 8, height = 6)
  ggsave(plot = traceplot1b, f1b, width = 8, height = 6)
  ggsave(plot = traceplot2b, f2b, width = 8, height = 6)
  
  
  traceplot0s <- plot_traces_ver(species, 7500, 7600, 'WT/WT')
  traceplot1s <- plot_traces_ver(species, 107500, 107600, 'WT/MUT')
  traceplot2s <- plot_traces_ver(species, 207500, 207600, 'MUT/MUT')
  
  f0s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_short.pdf')
  f1s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_short.pdf')
  f2s<-paste0(plotdir, '../../v1.6.2and5/panel_drafts/1cd_example_traceAndHisto/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_short.pdf')
  
  ggsave(plot = traceplot0s, f0s, width = 8, height = 6)
  ggsave(plot = traceplot1s, f1s, width = 8, height = 6)
  ggsave(plot = traceplot2s, f2s, width = 8, height = 6)
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.pdf')
  f1sb<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.pdf')
  f2sb<-paste0(plotdir, 'traces/lowCV_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.pdf')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter(time %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir, 'traces/lowCV_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir, 'traces/lowCV_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.pdf'))
  
}


# change from wt/wt to wt/mut analysis

sumstats_wt_vs_het_all <- ggplot(allstats_full %>% 
         group_by(version, paramset, product, mutated_alleles) %>% 
         pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:HDTpval) %>% 
         pivot_wider(names_from = mutated_alleles, values_from = value) %>% filter(product == 'B1'), 
       aes(`0`, `1`)) + 
  geom_point(alpha = 0.5, stroke = 0) + 
  facet_wrap(~stat, scales = 'free') +
  theme_bw() +
  xlab('WT/WT') +
  ylab('WT/MUT') +
  ggtitle('summary stats in wt/wt and heterozygous genotypes')
ggsave(sumstats_wt_vs_het_all, file = paste0(plotdir, 'sumstats_wt_vs_het_all.pdf'))

hetmean10set <- compared_stats %>% filter(product == 'B1', compare == 'lfc21', mean_denom > 10) %>% dplyr::select(version, paramset) %>% unique()

sumstats_wt_vs_het_hetmean10 <- ggplot(allstats_full %>% 
                                   group_by(version, paramset, product, mutated_alleles) %>% 
                                   pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:HDTpval) %>% 
                                   pivot_wider(names_from = mutated_alleles, values_from = value) %>% filter(product == 'B1') %>%
                                     inner_join(hetmean10set, by=c('version', 'paramset')), 
                                 aes(`0`, `1`)) + 
  geom_point(alpha = 0.5, stroke = 0) + 
  facet_wrap(~stat, scales = 'free') +
  theme_bw() +
  xlab('WT/WT') +
  ylab('WT/MUT') +
  ggtitle('summary stats in wt/wt and heterozygous genotypes\nMinimum mean in het = 10')
ggsave(sumstats_wt_vs_het_hetmean10, file = paste0(plotdir, 'sumstats_wt_vs_het_hetmean10.pdf'))


unistats<-unique(compared_stats$stat)
# 
# for (st in unistats) {
#   
#   p1 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)))) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Change in mean')
#   ggsave(p1, file = paste0(plotdir, 'changein_', st, '_vs_mean.pdf'), width = 16, height = 8)
#   
#   p2 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_denom, eval(as.symbol(st)))) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Starting mean abundance')
#   
#   p3 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)), color = mean_denom)) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Change in mean')
#   
#   p1l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)))) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Change in mean')
#   ggsave(p1, file = paste0(plotdir, 'changein_', st, '_vs_mean.pdf'), width = 16, height = 8)
#   
#   p2l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_denom, eval(as.symbol(st)))) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Starting mean abundance')
#   
#   p3l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)), color = mean_denom)) + 
#     geom_point() + 
#     geom_text(aes(label = paramset)) +
#     facet_grid(product~compare, scales = 'free') +
#     theme_bw() +
#     ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
#     ylab(paste0('Change in ', st)) +
#     xlab('Change in mean')
#   
#   ggsave(p1, file = paste0(plotdir, 'Diff_changein_', st, '_vs_changeinmean.pdf'), width = 16, height = 8)
#   ggsave(p2, file = paste0(plotdir, 'Diff_changein_', st, '_vs_startingmean.pdf'), width = 16, height = 8)
#   ggsave(p3, file = paste0(plotdir, 'Diff_changein_', st, '_vs_changeinmean_colorstartingmean.pdf'), width = 16, height = 8)
#   
#   ggsave(p1l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_changeinmean.pdf'), width = 16, height = 8)
#   ggsave(p2l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_startingmean.pdf'), width = 16, height = 8)
#   ggsave(p3l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_changeinmean_colorstartingmean.pdf'), width = 16, height = 8)
#   
# }
# 
# lfc10_lhs <- compared_stats %>%
#   filter(compare == 'lfc10', product %in% c('A1', 'B1')) %>%
#   dplyr::select(product, paramset, stat, diff, mean_denom) %>%
#   pivot_wider(names_from = stat, values_from = diff) %>%
#   inner_join(lhs_sets, by = 'paramset') 
# pdf(paste0(plotdir, 'corrplot_1vs0mut.pdf'), width = 10, height = 10)
# cor_B1_paramratio_stats10 <- corrplot.mixed(cor(as.matrix(lfc10_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
# dev.off()
# 
# lfc20_lhs <- compared_stats %>%
#   filter(compare == 'lfc20', product %in% c('A1', 'B1')) %>%
#   dplyr::select(product, paramset, stat, diff, mean_denom) %>%
#   pivot_wider(names_from = stat, values_from = diff) %>%
#   inner_join(lhs_sets, by = 'paramset') 
# 
# pdf(paste0(plotdir, 'corrplot_2vs0mut.pdf'), width = 10, height = 10)
# cor_B1_paramratio_stats20 <- corrplot.mixed(cor(as.matrix(lfc20_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
# dev.off()
# 
# lfc21_lhs <- compared_stats %>%
#   filter(compare == 'lfc21', product %in% c('A1', 'B1')) %>%
#   dplyr::select(product, paramset, stat, diff, mean_denom) %>%
#   pivot_wider(names_from = stat, values_from = diff) %>%
#   inner_join(lhs_sets, by = 'paramset') 
# 
# pdf(paste0(plotdir, 'corrplot_2vs1mut.pdf'), width = 10, height = 10)
# cor_B1_paramratio_stats21 <- corrplot.mixed(cor(as.matrix(lfc21_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
# dev.off()

# filter and re-do analyses
# mean_denom > 10
# Hill n < 5
lfc10_lhsf <- compared_stats %>%
  filter(compare == 'lfc10', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, version, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets_full, by = c('version', 'paramset')) %>%
  filter(mean_denom > 10, Hill_coefficient_n < 5)

pdf(paste0(plotdir, 'corrplot_1vs0mut_filt.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats10f <- corrplot.mixed(cor(as.matrix(lfc10_lhsf %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-c(product, version)))), tl.pos = 'lt')  
dev.off()

lfc20_lhsf <- compared_stats %>%
  filter(compare == 'lfc20', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') %>%
  filter(mean_denom > 10, Hill_coefficient_n < 5)

pdf(paste0(plotdir, 'corrplot_2vs0mut_filt.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats20f <- corrplot.mixed(cor(as.matrix(lfc20_lhsf %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
dev.off()

lfc21_lhsf <- compared_stats %>%
  filter(compare == 'lfc21', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') %>%
  filter(mean_denom > 10, Hill_coefficient_n < 5)

pdf(paste0(plotdir, 'corrplot_2vs1mut_filt.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats21f <- corrplot.mixed(cor(as.matrix(lfc21_lhsf %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
dev.off()




