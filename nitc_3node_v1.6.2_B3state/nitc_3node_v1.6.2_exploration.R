
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
datadir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2/samples/'
plotdir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2/samples/exploratory_analysis/'
setwd(datadir)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}
paramsets <- 1:9900
lhs_sets <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets %<>%
  mutate(paramset = 1:nrow(lhs_sets))

allstats <- list()
allparams <- list()
for (paramset in paramsets){
  
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
    ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_paramset_', as.character(paramset), '.pdf'))
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

  if(is.null(dim(allstats))) {
    allstats <- spstats
  } else {
    allstats %<>% bind_rows(spstats)
  }

  if(is.null(dim(allparams))) {
    allparams <- params
  } else {
    allparams %<>% bind_rows(params)
  }
  
}
write.csv(allstats, file = paste0(plotdir, 'summarystats.csv'))
# compare summary stats
pseud = 0.01

compared_stats <- allstats %>% 
  group_by(paramset, product, mutated_alleles) %>% 
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
  inner_join(allstats %>% 
               dplyr::select(mutated_alleles, product, paramset, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))


unistats<-unique(compared_stats$stat)

for (st in unistats) {
  
  p1 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)))) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Change in mean')
  ggsave(p1, file = paste0(plotdir, 'changein_', st, '_vs_mean.pdf'), width = 16, height = 8)
  
  p2 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_denom, eval(as.symbol(st)))) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Starting mean abundance')
  
  p3 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('delta', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)), color = mean_denom)) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Change in mean')
  
  p1l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)))) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Change in mean')
  ggsave(p1, file = paste0(plotdir, 'changein_', st, '_vs_mean.pdf'), width = 16, height = 8)
  
  p2l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_denom, eval(as.symbol(st)))) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Starting mean abundance')
  
  p3l <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product', grepl('lfc', compare)) %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)), color = mean_denom)) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Change in mean')
  
  ggsave(p1, file = paste0(plotdir, 'Diff_changein_', st, '_vs_changeinmean.pdf'), width = 16, height = 8)
  ggsave(p2, file = paste0(plotdir, 'Diff_changein_', st, '_vs_startingmean.pdf'), width = 16, height = 8)
  ggsave(p3, file = paste0(plotdir, 'Diff_changein_', st, '_vs_changeinmean_colorstartingmean.pdf'), width = 16, height = 8)
  
  ggsave(p1l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_changeinmean.pdf'), width = 16, height = 8)
  ggsave(p2l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_startingmean.pdf'), width = 16, height = 8)
  ggsave(p3l, file = paste0(plotdir, 'LFC_changein_', st, '_vs_changeinmean_colorstartingmean.pdf'), width = 16, height = 8)
  
}

lfc10_lhs <- compared_stats %>%
  filter(compare == 'lfc10', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') 
pdf(paste0(plotdir, 'corrplot_1vs0mut.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats10 <- corrplot.mixed(cor(as.matrix(lfc10_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
dev.off()

lfc20_lhs <- compared_stats %>%
  filter(compare == 'lfc20', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') 

pdf(paste0(plotdir, 'corrplot_2vs0mut.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats20 <- corrplot.mixed(cor(as.matrix(lfc20_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
dev.off()

lfc21_lhs <- compared_stats %>%
  filter(compare == 'lfc21', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') 

pdf(paste0(plotdir, 'corrplot_2vs1mut.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats21 <- corrplot.mixed(cor(as.matrix(lfc21_lhs %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
dev.off()

# filter and re-do analyses
# mean_denom > 10
# Hill n < 5
lfc10_lhsf <- compared_stats %>%
  filter(compare == 'lfc10', product %in% c('A1', 'B1')) %>%
  dplyr::select(product, paramset, stat, diff, mean_denom) %>%
  pivot_wider(names_from = stat, values_from = diff) %>%
  inner_join(lhs_sets, by = 'paramset') %>%
  filter(mean_denom > 10, Hill_coefficient_n < 5)

pdf(paste0(plotdir, 'corrplot_1vs0mut_filt.pdf'), width = 10, height = 10)
cor_B1_paramratio_stats10f <- corrplot.mixed(cor(as.matrix(lfc10_lhsf %>% filter(product == 'B1') %>% ungroup() %>% dplyr::select(-product))), tl.pos = 'lt')  
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



for (st in unistats) {
  
  pvs1 <- ggplot(allstats %>% inner_join(lhs_sets, by = 'paramset') %>% filter(product == 'B1', Hill_coefficient_n < 5)) + 
    geom_point(aes(mean_product, eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = 10), linetype = 2) +
    # geom_text(aes(mean_product, eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic() +
    ggtitle(paste0(st, ' vs mean, only including Hill n < 5'))
  
  pvs2 <- ggplot(allstats %>% inner_join(lhs_sets, by = 'paramset') %>% filter(product == 'B1', Hill_coefficient_n < 5)) + 
    geom_point(aes(log(mean_product), eval(as.symbol(st)))) +
    geom_vline(aes(xintercept = log(10)), linetype = 2) +
    # geom_text(aes(log(mean_product), eval(as.symbol(st)), label = as.character(paramset))) +
    facet_grid(~mutated_alleles) +
    theme_classic()  +
    ggtitle(paste0(st, ' vs mean, only including Hill n < 5'))
  
  ggsave(pvs1, file = paste0(plotdir, 'PerGenotype_', st, '_vs_mean.pdf'), width = 16, height = 8) 
  ggsave(pvs2, file = paste0(plotdir, 'PerGenotype_', st, '_vs_logmean.pdf'), width = 16, height = 8) 
}

# pick some traces: high and low bimodality, high and low CV
sampledir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2/samples/'
tracedir <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.2/fullTraces/'

high_lfc10_bimod <- lfc10_lhsf %>% filter(product == 'B1') %>% arrange(-bimodality_coef) %>% head(5)
low_lfc10_bimod <-  lfc10_lhsf %>% filter(product == 'B1', bimodality_coef > -0.01) %>% arrange(bimodality_coef) %>% head(5)
high_lfc10_cv <- lfc10_lhsf %>% filter(product == 'B1') %>% arrange(-cv_product) %>% head(5)
low_lfc10_cv <- lfc10_lhsf %>% filter(product == 'B1', cv_product > -0.01) %>% arrange(cv_product) %>% head(5)

for (pset in high_lfc10_bimod$paramset) {
  
  species<-as_tibble(read.csv(paste0(tracedir, 'initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/increasedBimod_0to1_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/increasedBimod_0to1_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/increasedBimod_0to1_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  species_sample <- species %>%
    mutate(paramset = pset) %>%
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
    ggsave(dist_plot, file = paste0(plotdir, 'traces/increasedBimod_0to1_distributions_q300_paramset_', as.character(pset), '.pdf'))
  
}


for (pset in low_lfc10_bimod$paramset) {
  
  species<-as_tibble(read.csv(paste0(tracedir, 'initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/stableBimod_0to1_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/stableBimod_0to1_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/stableBimod_0to1_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  species_sample <- species %>%
    mutate(paramset = pset) %>%
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
  ggsave(dist_plot, file = paste0(plotdir, 'traces/stableBimod_0to1_distributions_q300_paramset_', as.character(pset), '.pdf'))
  
}


for (pset in high_lfc10_cv$paramset) {
  
  species<-as_tibble(read.csv(paste0(tracedir, 'initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/increasedCV_0to1_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/increasedCV_0to1_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/increasedCV_0to1_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  species_sample <- species %>%
    mutate(paramset = pset) %>%
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
  ggsave(dist_plot, file = paste0(plotdir, 'traces/increasedCV_0to1_distributions_q300_paramset_', as.character(pset), '.pdf'))
  
}


for (pset in low_lfc10_cv$paramset) {
  
  species<-as_tibble(read.csv(paste0(tracedir, 'initialsim_species',as.character(pset),'.csv'), header = T)) %>%
    mutate(time = 1:300000)
  
  traceplot0 <- plot_traces(species, 7500, 8000, 'WT/WT')
  traceplot1 <- plot_traces(species, 107500, 108000, 'WT/MUT')
  traceplot2 <- plot_traces(species, 207500, 208000, 'MUT/MUT')
  
  f0<-paste0(plotdir, 'traces/stableCV_0to1_paramset', as.character(pset),'_WTWT.pdf')
  f1<-paste0(plotdir, 'traces/stableCV_0to1_paramset', as.character(pset),'_WTMUT.pdf')
  f2<-paste0(plotdir, 'traces/stableCV_0to1_paramset', as.character(pset),'_MUTMUT.pdf')
  
  ggsave(plot = traceplot0, f0, width = 8, height = 6)
  ggsave(plot = traceplot1, f1, width = 8, height = 6)
  ggsave(plot = traceplot2, f2, width = 8, height = 6)
  
  species_sample <- species %>%
    mutate(paramset = pset) %>%
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
  ggsave(dist_plot, file = paste0(plotdir, 'traces/stableCV_0to1_distributions_q300_paramset_', as.character(pset), '.pdf'))
  
}
