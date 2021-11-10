library(tidyverse)
library(magrittr)


setwd('~/code/grn_nitc/nitc_3node_v1.1/')

# load parameter sets

# load and process data
n_sets <- 100
for (i in 1:n_sets){
results <- t(as.matrix(read.csv('S_outpar_test_ifs.csv', header=F)))
colnames(results) <- c('orig_1', 'para_1', 'targ_1', 'targ_off', 'targ_on', 'para_off', 'para_on', 'orig_on', 'orig_off')
results %<>% as.data.frame() %>% as_tibble()
results[,'time'] <- 1:nrow(results)
results[,'paramset'] <- i

# sample every 300 timesteps for pseudo-single-cells pre-mutation (0-10k - "pre") and post-mutation (10k-20k "post imm"; 20k-30k "post lat")


# calculate gini


# store

}
orig_color = 'black'
para_color = 'gray'
targ_color = 'blue'
nt = 300

## Load 