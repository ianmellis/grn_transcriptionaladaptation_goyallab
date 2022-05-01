library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(diptest)
library(e1071)
library(ggrepel)
library(corrplot)
library(tibble)
library(svglite)
library()
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

lhs_sets2_Hn5<- lhs_sets2 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets2 <- lhs_sets2_Hn5$paramset

# calculate stats
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
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
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
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  
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
  


# temp: collate all data
setwd(datadir2)
all_species_q300 <- list()
for (paramset in paramsets2){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.2') %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, version, time, mutated_alleles)
  
  
  if(is.null(dim(all_species_q300))) {
    all_species_q300 <- species_sample
  } else {
    all_species_q300 %<>% bind_rows(species_sample)
  }
  
}

setwd(datadir5)
for (paramset in paramsets5){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.5') %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, version, time, mutated_alleles)
  
  
  if(is.null(dim(all_species_q300))) {
    all_species_q300 <- species_sample
  } else {
    all_species_q300 %<>% bind_rows(species_sample)
  }
  
}
write.csv(all_species_q300, file = paste0('/Volumes/IAMYG1/grn_nitc_data/v1.6.2and5/all_species_q300.csv'))


# summary stats on all data from the collated table, without reloading each file
allstats_full1 <- list()
psets <- data.frame(
  version = as.character(c(rep('1.6.2', length(paramsets2)), rep('1.6.5', length(paramsets5)))),
  paramset = c(paramsets2, paramsets5)
)
for (i in 1:nrow(psets)){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  ver = as.character(psets[i,'version'])
  paramset1 = psets[i,'paramset']
  
  species <- all_species_q300 %>% ungroup() %>%
    filter(version == ver & paramset == paramset1)
  
  species_sample <- species %>%
    pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')
  
  paramset = paramset1
  if(paramset %% 100 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_v1.6.5_paramset_', as.character(paramset), '_2.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  
  if(is.null(dim(allstats_full1))) {
    allstats_full1 <- spstats
  } else {
    allstats_full1 %<>% bind_rows(spstats)
  }
  
}



# pull specific parameter sets and plot histograms and traces for figure (svg files)

setsToPlot <- tibble(
  version = c('1.6.2', '1.6.5', '1.6.2', '1.6.5', '1.6.5', '1.6.2', '1.6.2', '1.6.2', '1.6.2', 
              '1.6.5', '1.6.5', '1.6.5', '1.6.5', '1.6.5', '1.6.5', '1.6.2', '1.6.2', '1.6.2', '1.6.2'),
  paramset = c(7000, 4900, 9400, 3000, 8100, 2800, 3300, 5800, 9900, 600, 3500, 
               4365, 3664, 4738, 7118, 4274, 574, 5873, 4290),
  classID = c('exponential', 'bimodal', 'gaussian', 'uniform', 'heavyTail', 'uniform', 'uniform', 'uniform', 'uniform', 'uniform', 'uniform', 'gaussian', 'gaussian', 'gaussian', 'gaussian', 'uniform', 'uniform', 'uniform', 'uniform')
)

if(!dir.exists(paste0(plotdir, 'panel_drafts/'))){
  dir.create(paste0(plotdir, 'panel_drafts/'))
}

for (ind in 1:nrow(setsToPlot)){
  
  pset = setsToPlot$paramset[ind]
  ver = setsToPlot$version[ind]
  classlab = setsToPlot$classID[ind]
  
  if(!dir.exists(paste0(plotdir, 'panel_drafts/', classlab))){
    dir.create(paste0(plotdir, 'panel_drafts/', classlab))
  }
  if(!dir.exists(paste0(plotdir, 'panel_drafts/', classlab, '/if_panel_1C/'))){
    dir.create(paste0(plotdir, 'panel_drafts/', classlab, '/if_panel_1C/'))
  }  
  if(!dir.exists(paste0(plotdir, 'panel_drafts/', classlab, '/supp/'))){
    dir.create(paste0(plotdir, 'panel_drafts/', classlab, '/supp/'))
  }
  if(ver == '1.6.2') {
    tracedir = datadir2
  } else {
    tracedir = datadir5
  }
  setwd(tracedir)
  species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0n <- plot_traces_ver_nonons(species, 7500, 8000, 'WT/WT')
  traceplot1n <- plot_traces_ver_nonons(species, 107500, 108000, 'WT/MUT')
  traceplot2n <- plot_traces_ver_nonons(species, 207500, 208000, 'MUT/MUT')
  
  f0n<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_noNonsense.svg')
  f1n<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_noNonsense.svg')
  f2n<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_noNonsense.svg')
  
  ggsave(plot = traceplot0n, f0n, width = 8, height = 6)
  ggsave(plot = traceplot1n, f1n, width = 8, height = 6)
  ggsave(plot = traceplot2n, f2n, width = 8, height = 6)
  
  
  traceplot0 <- plot_traces_ver(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces_ver(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces_ver(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'panel_drafts/', classlab, '/if_panel_1C/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT.svg')
  f1<-paste0(plotdir, 'panel_drafts/', classlab, '/if_panel_1C/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT.svg')
  f2<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT.svg')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  
  traceplot0sb <- plot_traces_ver_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1sb <- plot_traces_ver_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2sb <- plot_traces_ver_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0sb<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_Bfocus_short.svg')
  f1sb<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_Bfocus_short.svg')
  f2sb<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_Bfocus_short.svg')
  
  ggsave(plot = traceplot0sb, f0sb, width = 8, height = 6)
  ggsave(plot = traceplot1sb, f1sb, width = 8, height = 6)
  ggsave(plot = traceplot2sb, f2sb, width = 8, height = 6)
  
  traceplot0nb <- plot_traces_ver_nonons_Bfocus(species, 7500, 7600, 'WT/WT')
  traceplot1nb <- plot_traces_ver_nonons_Bfocus(species, 107500, 107600, 'WT/MUT')
  traceplot2nb <- plot_traces_ver_nonons_Bfocus(species, 207500, 207600, 'MUT/MUT')
  
  f0nb<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTWT_noNonsense_Bfocus_short.svg')
  f1nb<-paste0(plotdir, 'panel_drafts/', classlab, '/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_WTMUT_noNonsense_Bfocus_short.svg')
  f2nb<-paste0(plotdir, 'panel_drafts/', classlab, '/supp/trace_in_wtmut_version_',as.character(ver), '_paramset', as.character(pset),'_MUTMUT_noNonsense_Bfocus_short.svg')
  
  ggsave(plot = traceplot0nb, f0nb, width = 8, height = 6)
  ggsave(plot = traceplot1nb, f1nb, width = 8, height = 6)
  ggsave(plot = traceplot2nb, f2nb, width = 8, height = 6)
  
  
  
  species_sample <- species %>%
    mutate(paramset = pset, version = ver) %>%
    filter( (time+1) %% 300 == 0, (time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
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
    geom_rug() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot, file = paste0(plotdir, 'panel_drafts/', classlab, '/supp/histogram_in_wtmut_distributions_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
  dist_plot2<-ggplot(species_sample %>% filter(mutated_alleles < 2), aes(abundance)) +
    geom_histogram() +
    geom_rug() +
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot2, file = paste0(plotdir, 'panel_drafts/', classlab, '/if_panel_1C/histogram_in_wtmut_distributions_WTWT_WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
  dist_plot3 <- ggplot(species_sample %>% filter(mutated_alleles == 1, product == 'B1'), aes(abundance)) +
    geom_histogram() +
    geom_rug() +
    ggtitle(paste0('Version ', as.character(ver), ' Parameter set ', as.character(pset))) +
    theme_classic()
  ggsave(dist_plot3, file = paste0(plotdir, 'panel_drafts/', classlab, '/histogram_in_wtmut_distribution_B1WTMUTonly_q300_version_',as.character(ver), '_paramset_', as.character(pset), '.svg'))
  
}

# summary stats draft1
allstats_full <- bind_rows(allstats2 %>% mutate(version = '1.6.2'), allstats5 %>% mutate(version = '1.6.5'))
allparams_full <- bind_rows(allparams2 %>% mutate(paramset = 1:9900, version = '1.6.2'),
                            allparams5 %>% mutate(paramset = 1:10000, version = '1.6.5'))
lhs_sets_full <- bind_rows(lhs_sets2 %>% mutate(version = '1.6.2'), lhs_sets5 %>% mutate(version = '1.6.5'))

pseud = 0.01

compared_stats <- allstats_full %>% 
  group_by(version, paramset, product, mutated_alleles) %>% 
  pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:entropy90) %>% 
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


## Normalize stats to mean

loess_fitted_allstats_2 <- allstats2 %>% mutate(version = '1.6.2')
for (stat in unistats[unistats != 'mean_product']) {
  
  statdat <- list()
  
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    # for (ma in 0:2) {
      
      tempdat <- allstats2 %>% mutate(version = '1.6.2') %>% # remove for full
        filter(#mutated_alleles == ma,
               product == gene) %>%
        dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat)
      
      loess1 <- loess(eval(as.symbol(stat)) ~ mean_product, data = tempdat, span = 0.1)
      
      l1dat <- data.frame(mean_product = loess1$x,
                          stat = loess1$fitted,
                          resid = loess1$residuals,
                          version = tempdat$version,
                          paramset = tempdat$paramset,
                          product = gene,
                          mutated_alleles = tempdat$mutated_alleles)
      colnames(l1dat)[2] <- paste0(stat,'_fitted')
      colnames(l1dat)[3] <- paste0(stat,'_residual')
      
      lplot1 <- ggplot() +
        geom_point(data = tempdat, aes(mean_product, eval(as.symbol(stat))), alpha = 0.1) +
        geom_point(data = l1dat, aes(mean_product, eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
        theme_classic() +
        ylab(stat) +
        xlab('Mean') +
        ggtitle(paste0(stat, ' vs mean, with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
      
      lplot2 <- ggplot() +
        geom_point(data = tempdat, aes(log(mean_product), eval(as.symbol(stat))), alpha = 0.1) +
        geom_point(data = l1dat, aes(log(mean_product), eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
        theme_classic() +
        ylab(stat) +
        xlab('Log(Mean)') +
        ggtitle(paste0(stat, ' vs log(mean), with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
      
      ggsave(lplot1, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2only.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_v1.6.2only.pdf'), width = 5, height = 5)
      ggsave(lplot2, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2only.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_log_v1.6.2only.pdf'), width = 5, height = 5)
      
      if(is.null(dim(statdat))){
        statdat <- l1dat
      } else {
        statdat %<>% bind_rows(l1dat)
      }
    #}
    
  }
  
  statdat$version <- as.character(statdat$version)
  statdat$product <- as.character(statdat$product)
  statdat$paramset <- as.numeric(statdat$paramset)
  
  loess_fitted_allstats_2 %<>% left_join(as_tibble(statdat) %>% dplyr::select(-mean_product), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  statdat2 <- sliding_window_normalize(as_tibble(statdat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 25)
  
  loess_fitted_allstats_2 %<>% left_join(statdat2 %>% dplyr::select(-c('mean_product', paste0(stat,'_residual'), paste0(stat,'_fitted'))), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  td1<-allstats2 %>% mutate(version = '1.6.2') %>% # remove for full
    dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat)
  
  lplot_all_stat <- ggplot() +
    geom_point(data = td1 %>% filter(mean_product>10), aes(log(mean_product), eval(as.symbol(stat))), stroke=0, alpha = 0.05) +
    geom_density2d(data = td1 %>% filter(mutated_alleles == 0, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'blue') +
    geom_density2d(data = td1 %>% filter(mutated_alleles == 1, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'red') +
    geom_density2d(data = td1 %>% filter(mutated_alleles == 2, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'green') +
    theme_classic()
  
  lplot_all_statLOESS <-  ggplot() +
    geom_point(data = statdat %>% filter(mean_product>10), aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual')))), stroke=0, alpha = 0.05) +
    geom_density2d(data = statdat %>% filter(mutated_alleles == 0, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'blue') +
    geom_density2d(data = statdat %>% filter(mutated_alleles == 1, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'red') +
    geom_density2d(data = statdat %>% filter(mutated_alleles == 2, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'green') +
    theme_classic()
   
  lplot_all_statLOESSSWN <-  ggplot() +
    geom_point(data = statdat2 %>% filter(mean_product>10), aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual_swn')))), stroke=0, alpha = 0.05) +
    geom_density2d(data = statdat2 %>% filter(mutated_alleles == 0, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'blue') +
    geom_density2d(data = statdat2 %>% filter(mutated_alleles == 1, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'red') +
    geom_density2d(data = statdat2 %>% filter(mutated_alleles == 2, product == 'B1',mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'green') +
    theme_classic()
  
  pdf(paste0(paste0(plotdir, 'LOESSplots_', stat, 'vsLogMean_B1_v1.6.2only.pdf')), width = 10, height = 7)
  grid.arrange(lplot_all_stat,lplot_all_statLOESS,lplot_all_statLOESSSWN, ncol=3)
  dev.off()
  
}



bimfilt <- 0.15
high_bimodality <- loess_fitted_allstats_full %>% 
  filter(bimodality_residual > bimfilt)

unimodal_symmetric <- loess_fitted_allstats_full %>% 
  filter(bimodality_residual <= bimfilt) # and not exponential (skewness limit?)

unimodal_exponential <- loess_fitted_allstats_full %>% 
  filter(skewness > 1, skewness < 3) #skewness limit?

unimodal_subexponential <- loess_fitted_allstats_full %>% 
  filter(skewness >= 3) #skewness limit higher? Techinically skewness of exp = 2

entfilt <- 0.15
high_entropy <- loess_fitted_allstats_full %>% 
  filter(entropy_residual > entfilt)



