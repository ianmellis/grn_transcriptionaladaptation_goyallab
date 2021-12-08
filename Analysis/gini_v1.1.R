library(tidyverse)
library(magrittr)
library(ineq)

setwd('~/code/grn_nitc/nitc_3node_v1.1/')

n_sets <- 100

# load parameter sets
params <- as.data.frame(read.csv('R_outpar_NITC_1_1_1.csv', header = F))
colnames(params) <- c('A_prod1', 'Anonsense_prod1','Aprime_prod1','B_prod1',
                      'A_deg1','Anonsense_deg1','Aprime_deg1','B_deg1',
                      'B_ondep1','B_ondep_prime','Aprimenitc1','A_off1','Aprime_off1','B_off1',
                      'A_proddiff1','Aprime_proddiff1','B_proddiff1','onbasal_a1','onbasal_aprime1',
                      'onbasal_b1','kA1','kA1_nonsense','kAprime1','kB1',
                      'nA1','nA1_nonsense','nAprime1','nB1')
params[,'paramset'] <- 1:100

# load and process data
samp_stats <- list()
transition_samples <- list()
for (i in 1:n_sets){
  
  cat(paste0('Working on parameter set ', as.character(i), '...\n'))
  
  results <- t(as.matrix(read.csv(paste0('S_outpar_', as.character(i), '.csv'), header=F)))
  colnames(results) <- c('orig_1', 'nonsense_1', 'para_1', 'targ_1', 'targ_off', 'targ_on', 'para_off', 'para_on', 'orig_on', 'orig_off', 'is_mutated', 'not_mutated')
  results %<>% as.data.frame() %>% as_tibble()
  results[,'time'] <- 1:nrow(results)
  results[,'paramset'] <- i
  results %<>% filter(time <= 30000)
  
  transition_samp <- results %>% filter(time >= 9000, time <= 11000)
  if(is.null(dim(transition_samples))){
    transition_samples<-transition_samp
  } else {
    transition_samples %<>% bind_rows(transition_samp)
  }
  
  # sample every 300 timesteps for pseudo-single-cells pre-mutation (0-10k - "pre") and post-mutation (10k-20k "post imm"; 20k-30k "post lat")
  res_samps <- results %>%
    filter(time %% 300 == 0) %>%
    mutate(epoch = case_when(
      time < 10000 ~ 'pre',
      time > 10000 & time < 20000 ~ 'post imm',
      time > 20000 ~ 'post lat'
    ))
  
  # calculate CV and gini
  res_samps_stats <- res_samps %>%
    dplyr::select(epoch, orig_1, nonsense_1, para_1, targ_1) %>%
    pivot_longer(!epoch, names_to = 'gene', values_to = 'count') %>%
    group_by(epoch, gene) %>%
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

# plot CV and gini values for each parameter set
epoch_levels <- c('pre', 'post imm', 'post lat')
samp_stats$epoch1 <- factor(samp_stats$epoch, levels = epoch_levels)
samp_stats_tall <- samp_stats %>%
  pivot_longer(!c(epoch, epoch1, gene, paramset), names_to = 'statistic', values_to = 'value') %>%
  as_tibble()

stat_plot1 <- ggplot() +
  facet_grid(gene ~ statistic, scales = 'free_y') +
  geom_line(data = samp_stats_tall %>% filter(statistic != 'mean_count1'), aes(epoch1, value, group = paramset)) +
  theme_classic()

ggsave(stat_plot1, file = 'first100paramsets_CVgini.pdf', width = 7, height = 12)

stat_plot2 <- ggplot() +
  facet_grid(gene ~ statistic, scales = 'free_y') +
  geom_line(data = samp_stats_tall %>% filter(statistic == 'mean_count1'), aes(epoch1, log(value), group = paramset)) +
  theme_classic()

ggsave(stat_plot2, file = 'first100paramsets_logmean.pdf', width = 4, height = 12)



# isolate any parameter sets with weird values
# look for parameter trends in sets with increasing CV/gini
delta_stats <- samp_stats_tall %>%
  dplyr::select(-epoch1) %>%
  pivot_wider(names_from = c(epoch, statistic), values_from = value) %>%
  mutate(imm_delta_CV = `post imm_CV` - pre_CV,
         lat_delta_CV = `post lat_CV` - pre_CV,
         imm_delta_gini = `post imm_gini` - pre_gini,
         lat_delta_gini = `post lat_gini` - pre_gini)

delta_stats_sum <- delta_stats %>%
  dplyr::select(gene, paramset, imm_delta_CV, lat_delta_CV, imm_delta_gini, lat_delta_gini) %>%
  pivot_longer(!c(gene, paramset), names_to = 'delta_stat', values_to = 'value') %>%
  inner_join(params, by = 'paramset')

ggplot() +
  facet_wrap(.~gene, scales = 'free_y') +
  geom_boxplot(data = delta_stats_sum, aes(delta_stat, value))

ggplot(delta_stats_sum %>% filter(gene == 'targ_1', delta_stat == 'lat_delta_CV'),
       aes(B_ondep_prime, value)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('target dependency on paralog') +
  ylab('change in target CV after mutation')

ggplot(delta_stats_sum %>% filter(gene == 'targ_1', delta_stat == 'lat_delta_CV'),
       aes(B_ondep_prime, value)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab('target dependency on paralog') +
  ylab('change in target CV after mutation')

orig_color = 'black'
para_color = 'gray'
targ_color = 'blue'
nt = 300


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
  geom_point(data = samp_stats_wide %>% filter(gene == 'targ_1'), aes(deltalogMean_01, deltaCV_01))

deltaMeanVsdeltaginitarg_0_1mut <- ggplot()+
  theme_classic() + 
  geom_point(data = samp_stats_wide %>% filter(gene == 'targ_1'), aes(deltalogMean_01, gini_01))
