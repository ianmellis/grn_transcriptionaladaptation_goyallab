library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)

setwd('~/code/grn_nitc/nitc_TFpromoter_3node_v1.5/')

paramset = 5

params<-read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'),header = T)

species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'.csv'), header = T))%>%
  mutate(time = 1:nrow(species),
         paramset = paramset)
# colnames(species)<-c('A1',
#                      'Anonsense1',
#                      'Aprime1',
#                      'B1',
#                      'Promoter1_unbound_targ_allele1',
#                      'Promoter1_boundbyorig_targ_allele1',
#                      'Promoter1_boundbypara_targ_allele1',
#                      'Promoter1_unbound_targ_allele2',
#                      'Promoter1_boundbyorig_targ_allele2',
#                      'Promoter1_boundbypara_targ_allele2',
#                      'Burst1_on_targ_allele1',
#                      'Burst1_on_targ_allele2',
#                      'Burst1_on_para_allele1',
#                      'Burst1_on_para_allele2',
#                      'Burst1_on_orig_allele1',
#                      'Burst1_on_orig_allele2',
#                      'Burst1_is_mutated_allele1',
#                      'Burst1_is_mutated_allele2')
# species %<>% as_tibble() %>%
#   mutate(time = 1:nrow(species),
#          paramset = paramset)


transition_samp_1 <- species %>%
  filter(time > 97000, time < 103000)

transition_samp_2 <- species %>%
  filter(time > 197000, time < 203000)

steady_state_wtwt <- species %>%
  filter(time > 400, time < 100000) 


steady_state_sample_wtwt <- species %>%
  filter(time > 400, time < 100000, (time + 1) %% 100 == 0) 

steady_state_sample_wtmut <- species %>%
  filter(time > 100400, time < 200000, (time + 1) %% 100 == 0)

steady_state_sample_mutmut <- species %>%
  filter(time > 200400, (time + 1) %% 100 == 0)

autocor_ww<-acf(steady_state_wtwt$B1, lag.max = 1000)

autocor_wws<-acf(steady_state_sample_wtwt$B1)
autocor_wms<-acf(steady_state_sample_wtmut$B1)
autocor_mms<-acf(steady_state_sample_mutmut$B1)

grid.arrange(plot(autocor_wws), plot(autocor_wms), plot(autocor_mms))

tempacf <- tibble(
  acf = autocor$acf,
  lag = autocor$lag,
  paramset = i
)

orig_color = 'black'
nons_color = 'firebrick2'
para_color = 'gray'
targ_color = 'dodgerblue2'
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


spec_plot2 <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_2, aes(time, A1), color = orig_color) +
  geom_line(data = transition_samp_2, aes(time, Anonsense1), color = nons_color) +
  geom_line(data = transition_samp_2, aes(time, Aprime1), color = para_color) +
  geom_line(data = transition_samp_2, aes(time, B1), color = targ_color) +
  ylab('Abundance')+ theme_classic()

burstOn_plot2 <- ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_2 %>% filter(time < 200000), aes(time, Burst1_on_orig_allele1 + 0.1), color = nons_color) +
  geom_line(data = transition_samp_2 %>% filter(time < 200000), aes(time, Burst1_on_orig_allele2 + 0.09), color = orig_color) +
  geom_line(data = transition_samp_2 %>% filter(time >= 200000), aes(time, Burst1_on_orig_allele1 + 0.1), color = nons_color) +
  geom_line(data = transition_samp_2 %>% filter(time >= 200000), aes(time, Burst1_on_orig_allele2 + 0.09), color = nons_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_para_allele1 + 0.05), color = para_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_para_allele2 + 0.04), color = para_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_targ_allele1 + 0), color = targ_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_targ_allele2 - 0.01), color = targ_color) +
  geom_vline(data = transition_samp_2, aes(xintercept = 200000), color = nons_color, linetype = 2) +
  ylab('Burst status')+ theme_classic()

transition_samp_2_pro <- transition_samp_2 %>% 
  mutate(promoter1_state = case_when(
    Promoter1_unbound_targ_allele1 == 1 ~ 'unbound',
    Promoter1_boundbyorig_targ_allele1 == 1 ~ 'orig',
    Promoter1_boundbypara_targ_allele1 == 1 ~ 'para'
  ),
  promoter2_state = case_when(
    Promoter1_unbound_targ_allele2 == 1 ~ 'unbound',
    Promoter1_boundbyorig_targ_allele2 == 1 ~ 'orig',
    Promoter1_boundbypara_targ_allele2 == 1 ~ 'para'
  ))
promoter_plot2 <- ggplot() +
  geom_linerange(data = transition_samp_2_pro, stat = 'identity', aes(time, ymin = 0.7, ymax = 1.3, color = promoter1_state)) +
  geom_linerange(data = transition_samp_2_pro, stat = 'identity', aes(time, ymin = 1.7, ymax = 2.3, color = promoter2_state)) +
  scale_color_manual(values = c('orig' = orig_color, 'para' = para_color, 'unbound' = 'white')) + theme_classic() + theme(legend.position = 'none') +
  ylab('Target allele promoter')

t2_plot <- grid.arrange(spec_plot2, burstOn_plot2, promoter_plot2, ncol=1)



 ggplot() +
  theme_classic() +
  geom_line(data = transition_samp_2 %>% filter(time < 200000), aes(time, Burst1_on_orig_allele1 + 0.1), color = nons_color) +
  geom_line(data = transition_samp_2 %>% filter(time < 200000), aes(time, Burst1_on_orig_allele2 + 0.09), color = orig_color) +
  geom_line(data = transition_samp_2 %>% filter(time >= 200000), aes(time, Burst1_on_orig_allele1 + 0.1), color = nons_color) +
  geom_line(data = transition_samp_2 %>% filter(time >= 200000), aes(time, Burst1_on_orig_allele2 + 0.09), color = nons_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_para_allele1 + 0.05), color = para_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_para_allele2 + 0.04), color = para_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_targ_allele1 + 0), color = targ_color) +
  geom_line(data = transition_samp_2, aes(time, Burst1_on_targ_allele2 - 0.01), color = targ_color) +
  geom_vline(data = transition_samp_2, aes(xintercept = 200000), color = nons_color, linetype = 2) +
  ylab('Burst status') + xlim(c(199700, 200300))
 
 ggplot() +
   theme_classic() +
   geom_line(data = transition_samp_2, aes(time, A1), color = orig_color) +
   geom_line(data = transition_samp_2, aes(time, Anonsense1), color = nons_color) +
   geom_line(data = transition_samp_2, aes(time, Aprime1), color = para_color) +
   geom_line(data = transition_samp_2, aes(time, B1), color = targ_color) +
   ylab('Abundance') + xlim(c(199700, 200300))
