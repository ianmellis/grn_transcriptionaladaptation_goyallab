## Utility functions for analysis of GRN output

library(e1071)

# calculate_bimodality calculates Sarle's bimodality coefficient for a finite sample (https://en.wikipedia.org/wiki/Multimodal_distribution#Bimodality_coefficient).
# Note that a unimodal distribution will tend toward 0, a uniform distribution will give a value of ~0.555, and a Bernoulli distribution will give a value of 1.
# Accepts:
# - x: a numeric vector.
#
# Returns:
# - bc: a numeric value of Sarle's bimodality coefficient, between 0 and 1. Will return 0 for spike-and-slab dist (e.g., when all values 0)
calculate_bimodality <- function(x) {
  
  nx <- length(x)
  
  if(sum(x == mean(x)) == nx) {
    bc = 0
  } else {
    
    sx <- skewness(x)
    kx <- kurtosis(x) # excess kurtosis = m_4 / m_2^2 - 3
    
    bc <- (sx^2 + 1) / (kx + (3*(nx-1)^2)/((nx-2)*(nx-3)))
  }
  return(bc)
  
}

# plot_traces plots example traces for data from a simulation given timepoint bounds.
# Accepts:
# - dat: a tibble with columns: time, A1, Anonsense1, Aprime1, B1, paramset
# - start: numeric, starting time
# - end: numeric, ending time
# - genostring: string, describing genotype, e.g., 'WT/WT'
# 
# Returns:
# - p1: ggplot with traces only
plot_traces <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.5) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.5) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.5) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.5) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.5) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}

plot_traces_ver_Bfocus <- function(dat, start, end, genostring){
  
  orig_color = 'black'
  nons_color = 'firebrick2'
  para_color = 'gray'
  targ_color = 'dodgerblue2'
  
  species <- dat %>% filter(time >= start, time <= end)
  
  spec_plot <- ggplot() +
    theme_classic() +
    geom_line(data = species, aes(time, A1), color = orig_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Anonsense1), color = nons_color, alpha = 0.2) +
    geom_line(data = species, aes(time, Aprime1), color = para_color, alpha = 0.2) +
    geom_line(data = species, aes(time, B1), color = targ_color, alpha = 0.9) +
    ylab('Abundance') +
    ggtitle(paste0('Abundance over time, version ', as.character(ver), ', parameter set ', as.character(pset),'\nGenotype ', genostring, ' steady state'))
  
  return(spec_plot)
}
