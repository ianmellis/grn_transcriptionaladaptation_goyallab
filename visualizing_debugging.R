library(tidyverse)
library(magrittr)


setwd('~/code/grn_nitc/nitc_3node_v1/')

results <- t(as.matrix(read.csv('S_outpar_test_ifs.csv', header=F)))
colnames(results) <- c('orig_1', 'para_1', 'targ_1', 'targ_off', 'targ_on', 'para_off', 'para_on', 'orig_on', 'orig_off')
results %<>% as.data.frame() %>% as_tibble()
results[,'time'] <- 1:nrow(results)

orig_color = 'black'
para_color = 'gray'
targ_color = 'blue'
nt = 300
spec_plot <- ggplot() +
  theme_classic() +
  geom_line(data = results[1:nt,], aes(time, orig_1), color = orig_color) +
  geom_line(data = results[1:nt,], aes(time, para_1), color = para_color) +
  geom_line(data = results[1:nt,], aes(time, targ_1), color = targ_color) +
  ylab('Abundance')

burstOn_plot <- ggplot() +
  theme_classic() +
  geom_line(data = results[1:nt,], aes(time, orig_on + 0.05), color = orig_color) +
  geom_line(data = results[1:nt,], aes(time, para_on + 0.1), color = para_color) +
  geom_line(data = results[1:nt,], aes(time, targ_on - 0.1), color = targ_color) +
  ylab('Burst status')
grid.arrange(spec_plot, burstOn_plot, ncol=1)


burstOff_plot <- ggplot() +
  theme_classic() +
  geom_line(data = results, aes(time, orig_off), color = orig_color) +
  geom_line(data = results, aes(time, para_off), color = para_color) +
  geom_line(data = results, aes(time, targ_off), color = targ_color)



####

results <- t(as.matrix(read.csv('S_outpar_test_ifs_nonsenseSpecies.csv', header=F)))
colnames(results) <- c('orig_1', 'nonsense_1', 'para_1', 'targ_1', 'targ_off', 'targ_on', 'para_off', 'para_on', 'orig_on', 'orig_off', 'is_mutated', 'not_mutated')
results %<>% as.data.frame() %>% as_tibble()
results[,'time'] <- 1:nrow(results)

orig_color = 'black'
nons_color = 'orange'
para_color = 'gray'
targ_color = 'blue'
t_start = 10900
nt = 300
spec_plot <- ggplot() +
  theme_classic() +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, orig_1), color = orig_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, nonsense_1), color = nons_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, para_1), color = para_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, targ_1), color = targ_color) +
  ylab('Abundance')

burstOn_plot <- ggplot() +
  theme_classic() +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, orig_on + 0.1), color = orig_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, para_on + 0.05), color = para_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, targ_on + 0), color = targ_color) +
  geom_line(data = results[t_start:(t_start+nt),], aes(time, is_mutated - 0.05), color = nons_color) +
  ylab('Burst status')
grid.arrange(spec_plot, burstOn_plot, ncol=1)


