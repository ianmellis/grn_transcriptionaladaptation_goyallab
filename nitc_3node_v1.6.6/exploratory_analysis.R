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
datadir6 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.6/samples/'
plotdir6 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.6/exploratory_analysis/'
setwd(datadir6)

if(!dir.exists(plotdir6)){
  dir.create(plotdir6)
}
paramsets6 <- 1:9900
lhs_sets6 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets6 %<>%
  mutate(paramset = 1:nrow(lhs_sets6))

lhs_sets6_Hn5<- lhs_sets6 %>% filter(Hill_coefficient_n < 5)

# paramsets6 <- lhs_sets6_Hn5$paramset

# calculate stats
allstats6 <- list()
allparams6 <- list()
for (paramset in paramsets6){
  
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
  
  if(paramset %% 100 == 0 | paramset %in% c(1,2,3)) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir6, 'distributions_q300_v1.6.6_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              skewness = skewness(abundance + 0.01),
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
  
  if(is.null(dim(allstats6))) {
    allstats6 <- spstats
  } else {
    allstats6 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams6))) {
    allparams6 <- params
  } else {
    allparams6 %<>% bind_rows(params)
  }
  
}

all_species_q300 <- list()
for (paramset in paramsets6){
  
  if(paramset %% 100 == 0){cat(paste0('Working on ', as.character(paramset), '\n'))}
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.6') %>%
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
write.csv(all_species_q300, file = paste0('/Volumes/IAMYG1/grn_nitc_data/v1.6.6/all_species_q300.csv'))

pseud = 0.01

allstats6 %<>% mutate(skewness = ifelse(is.na(skewness), 0, skewness),
                      version = 'v1.6.6')

compared_stats <- allstats6 %>% 
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
  inner_join(allstats_full1 %>% 
               dplyr::select(mutated_alleles, product, paramset, version, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset', 'version')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))

# filter to Hill n < 5 - now already done for allstats_full1
#compared_stats %<>% ungroup() %>% inner_join(lhs_sets_fall %>% filter(Hill_coefficient_n < 5) %>% dplyr::select(version, paramset), by = c('version','paramset'))

# plot wt/mut summary stats against mean expression

unistats<-unique(compared_stats$stat)


loess_fitted_allstats_all <- allstats6 
for (stat in unistats[unistats != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  statdat <- list()
  statdat1 <- list()
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    # for (ma in 0:2) {
    
    tempdat <- allstats_full1 %>% 
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
    
    # ggsave(lplot1, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2and5.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_v1.6.2only.pdf'), width = 5, height = 5)
    # ggsave(lplot2, file = paste0(plotdir, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.2and5.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_log_v1.6.2only.pdf'), width = 5, height = 5)
    
    l1dat$version <- as.character(l1dat$version)
    l1dat$product <- as.character(l1dat$product)
    l1dat$paramset <- as.numeric(l1dat$paramset)
    
    cat('sliding window normalizing...\n') # do this per-gene over all genotypes, rather than all genes over all genotypes...
    l1dat1 <- sliding_window_normalize(as_tibble(l1dat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 50)
    
    if(is.null(dim(statdat))){
      statdat <- l1dat
      statdat1 <- l1dat1
    } else {
      statdat %<>% bind_rows(l1dat)
      statdat1 %<>% bind_rows(l1dat1)
    }
    
    
  }
  
  loess_fitted_allstats_all %<>% left_join(as_tibble(statdat) %>% dplyr::select(-mean_product), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  loess_fitted_allstats_all %<>% left_join(statdat1 %>% dplyr::select(-c('mean_product', paste0(stat,'_residual'), paste0(stat,'_fitted'))), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  cat(paste0('Done with ', stat, '\n'))
  
}

write.csv(loess_fitted_allstats_all, file = paste0(plotdir, 'loess_fitted_allstats_all_snwRadius100.csv'), quote = F, row.names = F)
