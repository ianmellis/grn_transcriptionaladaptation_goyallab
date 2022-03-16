

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
    theme_classic(); dist_plot
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

fullsim <- as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(paramset),'.csv'), header = T)) %>% mutate(time = 1:300000)
randtimes <- sample(-100:100, 332*3, replace = T)
indexinds <- c(seq(400, 99999, by = 300), seq(100400, 199999, by = 300), seq(200400, 299999, by = 300)) 
randinds<- randtimes + indexinds

species1 <- fullsim %>% filter(time %in% randinds, time < 100000)

species_sample1 <- fullsim %>%
  mutate(paramset = paramset) %>%
  filter(time %in% randinds) %>%
  mutate(mutated_alleles = case_when(
    time < 100001 ~ 0,
    time > 100000 & time < 200001 ~ 1,
    time > 200000 ~ 2
  )) %>%
  dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, time, mutated_alleles) %>%
  pivot_longer(cols = A1:B1, names_to = 'product', values_to = 'abundance')

dist_plot<-ggplot(species_sample1, aes(abundance)) +
  geom_histogram() +
  facet_grid(mutated_alleles~product) +
  ggtitle(paste0('Parameter set ', as.character(paramset))) +
  theme_classic(); dist_plot

if(!dir.exists('autocor/')){
  dir.create('autocor')
}
for (paramset in 1:100){
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  pdf(paste0('paramset', as.character(paramset), '_B1autocorr'))
  
}