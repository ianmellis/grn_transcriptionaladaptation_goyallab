library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')

datadir <- '/Volumes/IAMYG1/grn_nitc_data/v1.5.1/'
plotdir <- '~/code/grn_nitc/nitc_TFpromoter_3node_v1.5.1/exploratory_analysis/'
setwd(datadir)

if(!dir.exists(plotdir)){
  dir.create(plotdir)
}

paramsets <- 1:100
allstats <- list()
allparams <- list()
for (paramset in paramsets){
  
  cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T)) %>%
    mutate(paramset = paramset,
           ssA_wt = 2*r_onbasal_A1*r_prodon_A1/(r_deg_A1*(r_onbasal_A1+r_off_A1)),
           HBA_wt = r_bind_byA1_B1*(ssA_wt^n_A1)/(k_A1^n_A1 + ssA_wt^n_A1),
           boundA_wt = HBA_wt/(HBA_wt+r_unbind_byA1_B1),
           onB_wt = r_bound_byA1_B1*boundA_wt/(r_bound_byA1_B1*boundA_wt + r_off_B1),
           ssB_wt = 2*r_prodon_B1*onB_wt/r_deg_B1)
  
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
  
  dist_plot<-ggplot(species_sample, aes(abundance)) +
    geom_histogram() + 
    facet_grid(mutated_alleles~product) +
    ggtitle(paste0('Parameter set ', as.character(paramset))) +
    theme_classic()
  
  ggsave(dist_plot, file = paste0(plotdir, 'distributions_q300_paramset_', as.character(paramset), '.pdf'))
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              HDTpval = dip.test(abundance)$p.value)
  
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
compared_stats <- allstats %>% 
  group_by(paramset, product, mutated_alleles) %>% 
  pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:HDTpval) %>% 
  pivot_wider(names_from = mutated_alleles, values_from = value) %>% 
  mutate(lfc10 = log2(`1`/`0`), 
         delta10 = `1`-`0`,
         lfc21 = log2(`2`/`1`), 
         delta21 = `2`-`1`,
         lfc20 = log2(`2`/`0`), 
         delta20 = `2`-`0`) %>%
  dplyr::select(-c(`0`:`2`)) %>% 
  pivot_longer(names_to = 'compare', values_to = 'diff', cols = lfc10:delta20) 

unistats<-unique(compared_stats$stat)

for (st in unistats) {
  
  p1 <- ggplot(compared_stats %>% filter(stat == st | stat == 'mean_product') %>% pivot_wider(names_from = stat, values_from = diff), aes(mean_product, eval(as.symbol(st)))) + 
    geom_point() + 
    geom_text(aes(label = paramset)) +
    facet_grid(product~compare, scales = 'free') +
    theme_bw() +
    ggtitle(paste0('Change in ', st, ' vs. change in mean')) +
    ylab(paste0('Change in ', st)) +
    xlab('Change in mean')
  ggsave(p1, file = paste0(plotdir, 'changein_', st, '_vs_mean.pdf'), width = 16, height = 8)
  
}

ssA_plot <- ggplot(inner_join(allstats %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                filter(mutated_alleles == 0) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), allparams, by = 'paramset'), aes(A1, ssA_wt)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWild-type genotype, 100 parameter sets')

ssB_plot <- ggplot(inner_join(allstats %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                filter(mutated_alleles == 0) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), allparams, by = 'paramset'), aes(B1, ssB_wt)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWild-type genotype, 100 parameter sets')
pdf(paste0(plotdir, 'steadystate_wt_AB.pdf'), width = 10, height = 5)
ss_wt<-grid.arrange(ssA_plot, ssB_plot, ncol=2)
dev.off()


## het ode45 from matlab
steadystate_ode45_wtmut <- as_tibble(read.csv('../steady_state_ODE45_wtmut.csv', header = T))
ssA_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(A1, ss_A1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of A1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean A1') +
  ylab('ODE45 steady-state A1')

ssB_plot <- ggplot(inner_join(allstats %>% 
                                filter(mutated_alleles == 1) %>% 
                                dplyr::select(paramset, product, mean_product) %>% 
                                group_by(paramset) %>% 
                                pivot_wider(names_from = product, values_from = mean_product), steadystate_ode45_wtmut, by = 'paramset'), aes(B1, ss_B1)) +
  geom_point() +
  theme_bw() +
  ggtitle('Steady-state approximation of B1 vs simulation\nWT/MUT genotype, 100 parameter sets') +
  xlab('Simulated pseduo-single-cell mean B1') +
  ylab('ODE45 steady-state B1')

pdf(paste0(plotdir, 'steadystate_wtmut_AB.pdf'), width = 10, height = 5)
ss_wt<-grid.arrange(ssA_plot, ssB_plot, ncol=2)
dev.off()
