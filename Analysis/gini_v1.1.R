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
  pivot_longer(!c(epoch, epoch1, gene, paramset), names_to = 'statistic', values_to = 'value')

stat_plot1 <- ggplot() +
  facet_grid(gene ~ statistic, scales = 'free_y') +
  geom_line(data = samp_stats_tall %>% filter(statistic != 'mean_count1'), aes(epoch1, value, group = paramset)) +
  theme_classic()

ggsave(stat_plot1, file = 'first100paramsets_CVgini.pdf', width = 7, height = 12)

stat_plot2 <- ggplot() +
  facet_grid(gene ~ statistic, scales = 'free_y') +
  geom_line(data = samp_stats_tall %>% filter(statistic == 'mean_count1'), aes(epoch1, log(value), group = paramset)) +
  theme_classic()

ggsave(stat_plot1, file = 'first100paramsets_CVgini.pdf', width = 7, height = 12)


# isolate any parameter sets with weird values


orig_color = 'black'
para_color = 'gray'
targ_color = 'blue'
nt = 300

## Load 