library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)

setwd('~/code/grn_nitc/nitc_TFpromoter_3node_v1.5/')

params<-read.csv('initialsim_rates1.csv',header = T)

species<-t(as.matrix(read.csv('initialsim_species1.csv', header = F)))
colnames(species)<-c('A1',
                     'Anonsense1',
                     'Aprime1',
                     'B1',
                     'Promoter1_unbound_targ_allele1',
                     'Promoter1_boundbyorig_targ_allele1',
                     'Promoter1_boundbypara_targ_allele1',
                     'Promoter1_unbound_targ_allele2',
                     'Promoter1_boundbyorig_targ_allele2',
                     'Promoter1_boundbypara_targ_allele2',
                     'Burst1_on_targ_allele1',
                     'Burst1_on_targ_allele2',
                     'Burst1_on_para_allele1',
                     'Burst1_on_para_allele2',
                     'Burst1_on_orig_allele1',
                     'Burst1_on_orig_allele2',
                     'Burst1_is_mutated_allele1',
                     'Burst1_is_mutated_allele2')
species %<>% as_tibble() %>%
  mutate(time = 1:nrow(species),
         paramset = 1)


transition_samp_1 <- species %>%
  filter(time > 97000, time < 103000)

transition_samp_2 <- species %>%
  filter(time > 197000, time < 203000)


orig_color = 'black'
nons_color = 'orange'
para_color = 'gray'
targ_color = 'blue'
t_start = min(transition_samp_1$time)
t_end = max(transition_samp_1$time)
spec_plot <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_1, aes(time, A1), color = orig_color) +
  geom_line(data = transition_samp_1, aes(time, Anonsense1), color = nons_color) +
  geom_line(data = transition_samp_1, aes(time, Aprime1), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, B1), color = targ_color) +
  ylab('Abundance')

burstOn_plot <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_1 %>% filter(time < 100000), aes(time, Burst1_on_orig_allele1 + 0.1), color = orig_color) +
  geom_line(data = transition_samp_1 %>% filter(time < 100000), aes(time, Burst1_on_orig_allele2 + 0.09), color = orig_color) +
  geom_line(data = transition_samp_1 %>% filter(time >= 100000), aes(time, Burst1_on_orig_allele1 + 0.1), color = nons_color) +
  geom_line(data = transition_samp_1 %>% filter(time >= 100000), aes(time, Burst1_on_orig_allele2 + 0.09), color = orig_color) +
  geom_line(data = transition_samp_1, aes(time, Burst1_on_para_allele1 + 0.05), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, Burst1_on_para_allele2 + 0.04), color = para_color) +
  geom_line(data = transition_samp_1, aes(time, Burst1_on_targ_allele1 + 0), color = targ_color) +
  geom_line(data = transition_samp_1, aes(time, Burst1_on_targ_allele2 - 0.01), color = targ_color) +
  geom_vline(data = transition_samp_1, aes(xintercept = 100000), color = nons_color, linetype = 2) +
  ylab('Burst status')
grid.arrange(spec_plot, burstOn_plot, ncol=1)
