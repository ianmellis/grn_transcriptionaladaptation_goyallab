library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(grid)
library(diptest)
library(e1071)
library(ggrepel)
library(corrplot)
library(tibble)
library(svglite)
library(entropy)
library(ggalluvial)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')


# edit as needed
datadir51 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.1/samples/'
plotdir51 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.5.1/exploratory_analysis/'
setwd(datadir51)

if(!dir.exists(plotdir51)){
  dir.create(plotdir51)
}
paramsets51 <- 1:3
lhs_sets51 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets51 %<>%
  mutate(paramset = 1:nrow(lhs_sets51))

lhs_sets51_Hn5<- lhs_sets51 %>% filter(Hill_coefficient_n < 5, paramset <= 9900)

paramsets51 <- lhs_sets51_Hn5$paramset

# calculate stats
allstats51 <- list()
allparams51 <- list()
for (paramset in paramsets51){
  
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
  
  if(paramset %% 1 == 0) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir51, 'distributions_q300_v1.6.5.1_paramset_', as.character(paramset), '.pdf'))
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
  
  if(is.null(dim(allstats51))) {
    allstats51 <- spstats
  } else {
    allstats51 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams51))) {
    allparams51 <- params
  } else {
    allparams51 %<>% bind_rows(params)
  }
  
}
