library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)

setwd('~/code/grn_nitc/nitc_3node_v1.2/')

#### Diploid prototype
n_sets = 100
params <- as.data.frame(read.csv('R_outpar_NITC_dip_1_1_1.csv', header = F))
colnames(params) <- c('A_prod1', 'Anonsense_prod1','Aprime_prod1','B_prod1',
                      'A_deg1','Anonsense_deg1','Aprime_deg1','B_deg1',
                      'B_ondep1','B_ondep_prime','Aprimenitc1','A_off1','Aprime_off1','B_off1',
                      'A_proddiff1','Aprime_proddiff1','B_proddiff1','onbasal_a1','onbasal_aprime1',
                      'onbasal_b1','kA1','kA1_nonsense','kAprime1','kB1',
                      'nA1','nA1_nonsense','nAprime1','nB1')
params[,'paramset'] <- 1:100
samp_stats <- list()
transition_samples <- list()
for (i in 1:n_sets){
  
  cat(paste0('Working on parameter set ', as.character(i), '...\n'))
  
  results <- t(as.matrix(read.csv(paste0('S_outpar_dip_', as.character(i), '_q500.csv'), header=F)))
  colnames(results) <- c('orig_1', 'nonsense_1', 'para_1', 'targ_1', 'targ_allele1_off', 'targ_allele1_on', 'targ_allele2_off', 'targ_allele2_on', 
                         'para_allele1_off', 'para_allele1_on', 'para_allele2_off', 'para_allele2_on', 
                         'orig_allele1_on', 'orig_allele1_off', 'orig_allele2_on', 'orig_allele2_off', 
                         'allele1_is_mutated', 'allele1_not_mutated', 'allele2_is_mutated', 'allele2_not_mutated')
  results %<>% as.data.frame() %>% as_tibble()
  results[,'time'] <- seq(500,150000, by = 500)
  results[,'paramset'] <- i
  # results %<>% filter(time <= 30000)
  
  transition_results <- t(as.matrix(read.csv(paste0('S_outpar_dip_', as.character(i), '_transitions.csv'), header=F)))
  colnames(transition_results) <- c('orig_1', 'nonsense_1', 'para_1', 'targ_1', 'targ_allele1_off', 'targ_allele1_on', 'targ_allele2_off', 'targ_allele2_on', 
                                    'para_allele1_off', 'para_allele1_on', 'para_allele2_off', 'para_allele2_on', 
                                    'orig_allele1_on', 'orig_allele1_off', 'orig_allele2_on', 'orig_allele2_off', 
                                    'allele1_is_mutated', 'allele1_not_mutated', 'allele2_is_mutated', 'allele2_not_mutated')
  transition_results %<>% as.data.frame() %>% as_tibble()
  transition_results[,'time'] <- c(49000:51000,99000:101000)
  transition_results[,'paramset'] <- i
  transition_samp_1 <- transition_results %>% filter(time <= 51000)
  transition_samp_2 <- transition_results %>% filter(time >= 99000)
  if(is.null(dim(transition_samples))){
    transition_samples<-transition_samp
  } else {
    transition_samples %<>% bind_rows(transition_samp)
  }
  
  # sample every 500 timesteps for pseudo-single-cells
  res_samps <- results %>%
    mutate(mutated_alleles = case_when(
      time < 50001 ~ 0,
      time > 50001 & time < 100001 ~ 1,
      time > 100001 ~ 2
    ))
  
  # calculate CV and gini
  res_samps_stats <- res_samps %>%
    dplyr::select(mutated_alleles, orig_1, nonsense_1, para_1, targ_1) %>%
    pivot_longer(!mutated_alleles, names_to = 'gene', values_to = 'count') %>%
    group_by(mutated_alleles, gene) %>%
    summarise(mean_count1 = mean(count+1),
              CV = sd(count + 1)/mean(count + 1),
              gini = ineq(count + 1, type = 'Gini', na.rm = T)) %>%
    mutate(paramset = i)
  
  # store
  if(is.null(dim(samp_stats))){
    samp_stats <- res_samps_stats
  } else {
    samp_stats %<>% bind_rows(res_samps_stats)
  }
}

orig_color = 'black'
nons_color = 'orange'
para_color = 'gray'
targ_color = 'blue'
t_start = min(transition_samp_1$time)
t_end = max(transition_samp_1$time)
nt = 300
spec_plot <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_1, aes(time, orig_1), color = orig_color) +
  geom_line(data = transition_samp_1, aes(time, nonsense_1), color = nons_color) +
  geom_line(data = transition_samp_1, aes(time, para_1), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, targ_1), color = targ_color) +
  ylab('Abundance')

burstOn_plot <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_1 %>% filter(time < 50000), aes(time, orig_allele1_on + 0.1), color = orig_color) +
  geom_line(data = transition_samp_1 %>% filter(time < 50000), aes(time, orig_allele2_on + 0.09), color = orig_color) +
  geom_line(data = transition_samp_1 %>% filter(time >= 50000), aes(time, orig_allele1_on + 0.1), color = nons_color) +
  geom_line(data = transition_samp_1 %>% filter(time >= 50000), aes(time, orig_allele2_on + 0.09), color = orig_color) +
  geom_line(data = transition_samp_1, aes(time, para_allele1_on + 0.05), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, para_allele2_on + 0.04), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, targ_allele1_on + 0), color = targ_color) +
  geom_line(data = transition_samp_1, aes(time, targ_allele2_on - 0.01), color = targ_color) +
  geom_vline(data = transition_samp_1, aes(xintercept = 50000), color = nons_color, linetype = 2) +
  ylab('Burst status')
grid.arrange(spec_plot, burstOn_plot, ncol=1)


meanVsCV_targ <- ggplot()+
  theme_classic() + 
  geom_point(data = samp_stats %>% filter(gene == 'targ_1'), aes(log(mean_count1), CV)) +
  facet_wrap(~mutated_alleles)

samp_stats_wide <- samp_stats %>% 
  pivot_wider(names_from = mutated_alleles, values_from = c(mean_count1, CV, gini)) %>%
  mutate(deltalogMean_01 = log(mean_count1_1) - log(mean_count1_0),
         deltalogMean_02 = log(mean_count1_2) - log(mean_count1_1),
         deltalogMean_12 = log(mean_count1_2) - log(mean_count1_0),
         deltaCV_01 = CV_1 - CV_0,
         deltaCV_02 = CV_2 - CV_0,
         deltaCV_12 = CV_2 - CV_1,
         gini_01 = gini_1 - gini_0,
         gini_02 = gini_2 - gini_0,
         gini_12 = gini_2 - gini_1)

deltaMeanVsdeltaCVtarg_0_1mut <- ggplot()+
  theme_classic() + 
  geom_point(data = samp_stats_wide %>% inner_join(params) %>% filter(gene == 'targ_1', deltalogMean_01 > -0.8, B_ondep_prime/B_ondep1 < 2), aes(deltalogMean_01, deltaCV_01, color = log2(B_ondep_prime/B_ondep1))) +
  scale_color_continuous(type = "viridis")

deltaMeanVsdeltaginitarg_0_1mut <- ggplot()+
  theme_classic() + 
  geom_point(data = samp_stats_wide %>% filter(gene == 'targ_1'), aes(deltalogMean_01, gini_01))

deltaMeanVsNITC_0_1mut <- ggplot()+
  theme_classic() + 
  geom_point(data = samp_stats_wide %>% filter(gene == 'targ_1') %>% inner_join(params), aes(deltalogMean_01, Aprimenitc1))


## correlations between parameters and pseudo-single-cell gene expression
flattenCorrMatrix <- function(cormat, pmat) { # from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

total_results <- samp_stats_wide %>% inner_join(params)
cors_targ<-rcorr(as.matrix(total_results %>% filter(gene == 'targ_1') %>% dplyr::select(-gene)))
cors_targ_tall<-flattenCorrMatrix(cors_targ$r, cors_targ$P)

var_intr <- 'deltaCV_01'
cors_targ_tall %>% filter(row == var_intr | column == var_intr) %>% arrange(-abs(cor))
cors_targ_tall %>% filter(row == 'deltalogMean_01' | column == 'deltalogMean_01') %>% arrange(-abs(cor))

lm()
