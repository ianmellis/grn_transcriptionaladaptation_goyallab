library(tidyverse)
library(magrittr)
library(ineq)
library(Hmisc)
library(gridExtra)
library(grid)
library(diptest)
library(e1071)
library(ggrepel)
library(corrplot)
library(tibble)
library(svglite)
library(entropy)
library(ggalluvial)
source('~/code/grn_nitc/Functions/grn_analysis_utilities.R')


# edit as needed
datadir8 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.8/samples/'
plotdir8 <- '/Volumes/IAMYG1/grn_nitc_data/v1.6.8/exploratory_analysis/'
setwd(datadir8)

if(!dir.exists(plotdir8)){
  dir.create(plotdir8)
}
paramsets8 <- 1:9900
lhs_sets8 <- as_tibble(read.csv('latinhyp_sampledSets.csv')) 
lhs_sets8 %<>%
  mutate(paramset = 1:nrow(lhs_sets8))

lhs_sets8_Hn5<- lhs_sets8 %>% filter(Hill_coefficient_n < 5)

pset = 2
ver = '1.6.8'
species<-as_tibble(read.csv(paste0('../fullTraces/initialsim_species',as.character(pset),'.csv'), header = T)) %>%
  mutate(time = 1:300000)

traceplot0n <- plot_traces_ver_nonons(species, 7500, 8000, 'WT/WT')
traceplot1n <- plot_traces_ver_nonons(species, 107500, 108000, 'WT/MUT')
traceplot2n <- plot_traces_ver_nonons(species, 207500, 208000, 'MUT/MUT')

# calculate stats
allstats8 <- list()
allparams8 <- list()
for (paramset in paramsets8){
  
  # cat(paste0('Working on ', as.character(paramset), '\n'))
  
  params<-as_tibble(read.csv(paste0('initialsim_rates',as.character(paramset),'.csv'), header = T))
  
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
  
  if(paramset %% 100 == 0 | paramset %in% c(1,2,3)) {
    cat(paste0('Working on ', as.character(paramset), '\n'))
    dist_plot<-ggplot(species_sample, aes(abundance)) +
      geom_histogram() +
      facet_grid(mutated_alleles~product) +
      ggtitle(paste0('Parameter set ', as.character(paramset))) +
      theme_classic()
    ggsave(dist_plot, file = paste0(plotdir8, 'distributions_q300_v1.6.8_paramset_', as.character(paramset), '.pdf'))
  }
  
  spstats <- species_sample %>%
    group_by(mutated_alleles, product, paramset) %>%
    summarise(mean_product = mean(abundance),
              sd_product = sd(abundance),
              cv_product = sd_product/(mean_product + 0.01),
              fano_product = sd_product^2/(mean_product + 0.01),
              skewness = skewness(abundance + 0.01),
              gini_product = ineq(abundance + 1, type = 'Gini', na.rm = T),
              bimodality_coef = calculate_bimodality(abundance),
              entropy = entropy(discretize(abundance, numBins = 30, r=c(0,max(1,max(abundance))))), .groups = 'keep')
  
  entr_temp <- tibble(
    mutated_alleles = numeric(),
    product = character(),
    paramset = numeric(),
    entropy95 = numeric(),
    entropy90 = numeric()
  )
  
  for (ma in c(0,1,2)) {
    
    for (gene in c('A1', 'Aprime1', 'Anonsense1', 'B1')) {
      
      subs <- species_sample %>%
        filter(mutated_alleles == ma, product == gene)
      
      simdist <- subs$abundance
      
      # 
      # expFit <- fitdistr(simdist, 'exponential')
      # 
      # expKS <- ks.test(simdist, 'pexp', expFit$estimate)
      
      #filter to remove top and bottom 2.5% of values to assess entropy of distribution bulk
      
      nv <- length(simdist)
      simdistfilt95 <- simdist[order(simdist)[round(0.025*nv):round(0.975*nv)]]
      simdistfilt90 <- simdist[order(simdist)[round(0.05*nv):round(0.95*nv)]]
      
      trow <- tibble(
        mutated_alleles = ma,
        paramset = paramset,
        product = gene,
        entropy95 = entropy(discretize(simdistfilt95, numBins = 30, r=c(0,max(1,max(simdistfilt95))))),
        entropy90 = entropy(discretize(simdistfilt90, numBins = 30, r=c(0,max(1,max(simdistfilt90)))))
      )
      
      entr_temp %<>% bind_rows(trow)
      
    }
  }
  
  spstats %<>% inner_join(entr_temp, by = c('mutated_alleles', 'product', 'paramset'))
  
  if(is.null(dim(allstats8))) {
    allstats8 <- spstats
  } else {
    allstats8 %<>% bind_rows(spstats)
  }
  
  if(is.null(dim(allparams8))) {
    allparams8 <- params
  } else {
    allparams8 %<>% bind_rows(params)
  }
  
}
write.csv(allstats8, file = paste0(plotdir8, 'summary_stats_perGenotype.csv'))


all_species_q300 <- list()
for (paramset in paramsets8){
  
  if(paramset %% 100 == 0){cat(paste0('Working on ', as.character(paramset), '\n'))}
  
  species<-as_tibble(read.csv(paste0('initialsim_species',as.character(paramset),'_q300.csv'), header = T))
  
  species_sample <- species %>%
    mutate(paramset = paramset,
           version = '1.6.8') %>%
    filter((time > 400 & time < 100001) | (time > 100400 & time < 200001) | (time > 200400 & time < 300001)) %>%
    mutate(mutated_alleles = case_when(
      time < 100001 ~ 0,
      time > 100000 & time < 200001 ~ 1,
      time > 200000 ~ 2
    )) %>%
    dplyr::select(A1, Aprime1, Anonsense1, B1, paramset, version, time, mutated_alleles)
  
  
  if(is.null(dim(all_species_q300))) {
    all_species_q300 <- species_sample
  } else {
    all_species_q300 %<>% bind_rows(species_sample)
  }
  
}
write.csv(all_species_q300, file = paste0('/Volumes/IAMYG1/grn_nitc_data/v1.6.8/all_species_q300.csv'))

pseud = 0.01

allstats8 %<>% mutate(skewness = ifelse(is.na(skewness), 0, skewness),
                      version = 'v1.6.7')

compared_stats <- allstats8 %>% 
  group_by(version, paramset, product, mutated_alleles) %>% 
  pivot_longer(names_to = 'stat', values_to = 'value', cols = mean_product:entropy90) %>% 
  pivot_wider(names_from = mutated_alleles, values_from = value) %>% 
  mutate(lfc10 = log2((`1`+pseud)/(`0` + pseud)), 
         delta10 = `1`-`0`,
         lfc21 = log2((`2`+pseud)/(`1` + pseud)), 
         delta21 = `2`-`1`,
         lfc20 = log2((`2`+pseud)/(`0` + pseud)), 
         delta20 = `2`-`0`) %>%
  dplyr::select(-c(`0`:`2`)) %>% 
  pivot_longer(names_to = 'compare', values_to = 'diff', cols = lfc10:delta20) %>%
  inner_join(allstats6 %>% 
               dplyr::select(mutated_alleles, product, paramset, version, mean_product) %>%
               pivot_wider(names_from = mutated_alleles, values_from = mean_product), by = c('product', 'paramset', 'version')) %>%
  mutate(mean_denom = case_when(
    compare %in% c('lfc10', 'delta10', 'lfc20', 'delta20') ~ `0`,
    compare %in% c('lfc21', 'delta21') ~ `1`))

# filter to Hill n < 5 - now already done for allstats_full1
#compared_stats %<>% ungroup() %>% inner_join(lhs_sets_fall %>% filter(Hill_coefficient_n < 5) %>% dplyr::select(version, paramset), by = c('version','paramset'))

# plot wt/mut summary stats against mean expression

unistats<-unique(compared_stats$stat)


loess_fitted_allstats_all8 <- allstats7
for (stat in unistats[unistats != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  statdat <- list()
  statdat1 <- list()
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    # for (ma in 0:2) {
    
    tempdat <- allstats8 %>% 
      filter(#mutated_alleles == ma,
        product == gene) %>%
      dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat)
    
    loess1 <- loess(eval(as.symbol(stat)) ~ mean_product, data = tempdat, span = 0.1)
    
    l1dat <- data.frame(mean_product = loess1$x,
                        stat = loess1$fitted,
                        resid = loess1$residuals,
                        version = tempdat$version,
                        paramset = tempdat$paramset,
                        product = gene,
                        mutated_alleles = tempdat$mutated_alleles)
    colnames(l1dat)[2] <- paste0(stat,'_fitted')
    colnames(l1dat)[3] <- paste0(stat,'_residual')
    
    lplot1 <- ggplot() +
      geom_point(data = tempdat, aes(mean_product, eval(as.symbol(stat))), alpha = 0.1) +
      geom_point(data = l1dat, aes(mean_product, eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
      theme_classic() +
      ylab(stat) +
      xlab('Mean') +
      ggtitle(paste0(stat, ' vs mean, with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
    
    lplot2 <- ggplot() +
      geom_point(data = tempdat, aes(log(mean_product), eval(as.symbol(stat))), alpha = 0.1) +
      geom_point(data = l1dat, aes(log(mean_product), eval(as.symbol(paste0(stat,'_fitted')))), color = 'red') +
      theme_classic() +
      ylab(stat) +
      xlab('Log(Mean)') +
      ggtitle(paste0(stat, ' vs log(mean), with LOESS fit to mean\nGene product: ', gene))#, ', mutated alleles: ', as.character(ma)))
    
    ggsave(lplot1, file = paste0(plotdir8, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.8.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_v1.6.2only.pdf'), width = 5, height = 5)
    ggsave(lplot2, file = paste0(plotdir8, 'LOESS_', stat, 'vsMean_',gene,'_v1.6.8.pdf'), width = 5, height = 5)#'_mutAlleles',ma,'_log_v1.6.2only.pdf'), width = 5, height = 5)
    
    l1dat$version <- as.character(l1dat$version)
    l1dat$product <- as.character(l1dat$product)
    l1dat$paramset <- as.numeric(l1dat$paramset)
    
    cat('sliding window normalizing...\n') # do this per-gene over all genotypes, rather than all genes over all genotypes...
    l1dat1 <- sliding_window_normalize(as_tibble(l1dat) %>% filter(mean_product>10), 'mean_product', paste0(stat,'_residual'), 50)
    
    if(is.null(dim(statdat))){
      statdat <- l1dat
      statdat1 <- l1dat1
    } else {
      statdat %<>% bind_rows(l1dat)
      statdat1 %<>% bind_rows(l1dat1)
    }
    
    
  }
  
  loess_fitted_allstats_all7 %<>% left_join(as_tibble(statdat) %>% dplyr::select(-mean_product), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  loess_fitted_allstats_all7 %<>% left_join(statdat1 %>% dplyr::select(-c('mean_product', paste0(stat,'_residual'), paste0(stat,'_fitted'))), by = c('version', 'paramset', 'mutated_alleles', 'product'))
  
  cat(paste0('Done with ', stat, '\n'))
  
}

write.csv(loess_fitted_allstats_all8, file = paste0(plotdir8, 'loess_fitted_allstats_all_snwRadius100.csv'), quote = F, row.names = F)



for (stat in unistats[unistats != 'mean_product']) {
  
  cat(paste0('working on ', stat, '\n'))
  # statdat <- list()
  
  for (gene in c('A1', 'Anonsense1', 'Aprime1', 'B1')) {
    
    
    statdatA = loess_fitted_allstats_all8 %>%
      dplyr::select(mutated_alleles, product, version, paramset, mean_product, stat, as.symbol(paste0(stat, '_residual')), as.symbol(paste0(stat, '_residual_swn'))) %>%
      filter(product == gene, mean_product>10) 
    
    lplot_all_stat_con <- ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(stat))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2, mean_product>10), aes(log(mean_product),eval(as.symbol(stat))), color = 'green') +
      theme_classic()
    
    lplot_all_statLOESS_con <-  ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual')))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual')))), color = 'green') +
      theme_classic()
    
    lplot_all_statLOESSSWN_con <-  ggplot() +
      geom_point(data = statdatA, aes(log(mean_product), eval(as.symbol(paste0(stat,'_residual_swn')))), stroke=0, alpha = 0.05) +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 0,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'blue') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 1,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'red') +
      geom_density2d(data = statdatA %>% filter(mutated_alleles == 2,mean_product>10), aes(log(mean_product),eval(as.symbol(paste0(stat,'_residual_swn')))), color = 'green') +
      theme_classic()
    
    pdf(paste0(paste0(plotdir8, 'LOESSplots_', stat, 'vsLogMean_', gene,'_v1.6.8.pdf')), width = 10, height = 7)
    grid.arrange(lplot_all_stat_con,lplot_all_statLOESS_con,lplot_all_statLOESSSWN_con, ncol=3,
                 top = textGrob(paste0(stat, ' vs. log(mean_product), ', gene, '\nStat, LOESS residual, Squeezed LOESS residual (radius=50)'),gp=gpar(fontsize=20,font=3)))
    dev.off()
    
  }
}

anver <- 5 # increase minimum bimodality_residual filter and make left-skew filter more stringent
bimfilt <- 0.1
# high_bimodality <- loess_fitted_allstats_all %>% 
#   filter(bimodality_coef_residual > bimfilt)
# 
# unimodal_symmetric <- loess_fitted_allstats_all %>% 
#   filter(bimodality_coef_residual <= bimfilt, abs(skewness) < 1) # and not exponential (skewness limit?)
# 
# unimodal_exponential <- loess_fitted_allstats_all %>% 
#   filter(bimodality_coef_residual <= bimfilt, skewness > 1, skewness < 3) #skewness limit?
# 
# unimodal_subexponential <- loess_fitted_allstats_all %>% 
#   filter(bimodality_coef_residual <= bimfilt, skewness >= 3) #skewness limit higher? Techinically skewness of exp = 2

entfilt <- 0.15
# high_entropy <- loess_fitted_allstats_full %>% 
#   filter(entropy_residual > entfilt)

basic_class_assignment_all8 <- loess_fitted_allstats_all8 %>%
  mutate(class_assignment = case_when(
    mean_product < 10 ~ 'low-average',
    bimodality_coef_residual > bimfilt & bimodality_coef > 0.555 ~ 'bimodal',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & abs(skewness) < 1 ~ 'unimodal symmetric',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness >= 1 & skewness < 3 ~ 'exponential',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness >= 3 ~ 'subexponential',
    (bimodality_coef_residual <= bimfilt | bimodality_coef <= 0.555) & skewness <= -1 ~ 'left-skewed unimodal'
    
  )) 

basic_class_assignment_all_forSankey8 <- basic_class_assignment_all6 %>%
  dplyr::select(version, paramset, product, mutated_alleles, class_assignment) %>%
  group_by(version, paramset, product) %>%
  pivot_wider(names_from = mutated_alleles, values_from = class_assignment) %>%
  group_by(product, `0`, `1`, `2`) %>%
  summarise(Freq = length(`0`)) %>%
  ungroup() %>%
  mutate(alluvID = 1:length(product)) %>%
  pivot_longer(`0`:`2`, names_to = 'mutated_alleles', values_to = 'class_assignment')

# sample paramsets from the sankey flow to visually inspect accuracy of assignments/changes
set.seed(73245)
basic_class_assignment_all_forSankey_forsamples6 <- basic_class_assignment_all6 %>%
  dplyr::select(version, paramset, product, mutated_alleles, class_assignment)  

classes6 = unique(basic_class_assignment_all8$class_assignment)

if(!dir.exists(paste0(plotdir6, 'stats_class_assignment_check_classv', as.character(anver)))) {
  dir.create(paste0(plotdir6, 'stats_class_assignment_check_classv', as.character(anver)))
}

classes_sankey6 <- ggplot(basic_class_assignment_all_forSankey6, aes(x = mutated_alleles, y=Freq,
                                                                     stratum = class_assignment, alluvium = alluvID, fill = class_assignment, label = class_assignment)) +
  facet_grid(product~.) +
  geom_flow() +
  geom_stratum(alpha = 0.5) +
  geom_text(stat = 'stratum', size = 3) + 
  theme(legend.position = 'none')
ggsave(classes_sankey6, file = paste0(plotdir8, 'stats_class_assignment_check_classv', as.character(anver),'/classes_sankey.pdf'))

basic_class_assignment_all_forpie6 <- basic_class_assignment_all6 %>%
  group_by(mutated_alleles, product, class_assignment) %>%
  summarise(nSets = length(product))

classes_pies6 <- ggplot(basic_class_assignment_all_forpie6, aes(x="", y=nSets, fill=class_assignment)) +
  geom_bar(stat='identity', width=1, color='white') +
  coord_polar('y', start=0) +
  facet_grid(product~mutated_alleles) +
  theme_void() +
  ggtitle('classes of all distributions in v1.6.6\nparameter sets with Hill n < 5\neach gene in each genotype')
ggsave(classes_pies6, file = paste0(plotdir8, 'stats_class_assignment_check_classv', as.character(anver),'/classes_pies.pdf'))

# decision tree for classes
# separate by gene and genotype
# use parameters as features, use class ID as label
# then: don't separate by genotype, use parameters as feature and vector of class IDs as label
# then: don't separate by gene (at least jointly analyze A1 and B1), use parameters as features and vector of classIDs for A1 and B1 as labels
# consider trees in base R https://www.datacamp.com/tutorial/decision-trees-R and tidymodels https://bcullen.rbind.io/post/2020-06-02-tidymodels-decision-tree-learning-in-r/

# get data into format: each row is a paramset, with cols: sampled LHS (numeric) and classes (factor) of each gene in each genotype

basic_class_assignment_all8_wide <- basic_class_assignment_all8 %>% 
  dplyr::select(version, paramset, mutated_alleles, product, class_assignment) %>%
  pivot_wider(names_from = c('mutated_alleles', 'product'), values_from = 'class_assignment')

classes_for_trees <- inner_join(basic_class_assignment_all8_wide, lhs_sets8, by = 'paramset')

temp_class_for_tree <- classes_for_trees %>%
  dplyr::select(colnames(lhs_sets6), `1_B1`) %>% ungroup() %>%
  dplyr::select(-paramset) %>%
  mutate(is_bimodal = ifelse(`1_B1` == 'bimodal', T, F)) %>%
  dplyr::select(-`1_B1`)

temp_class_for_tree$is_bimodal <- as.factor(temp_class_for_tree$is_bimodal)

temp.tree <- ctree(is_bimodal ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + A1_Aprime_prodon_ratio + r_addon_byA1_B1 + r_onbasal_A1 + r_onbasal_Aprime_ratio, data = temp_class_for_tree)


# decision tree for B1 unimodal symm to unimodal symm vs unimodal symm to other


temp_class_for_tree <- classes_for_trees %>%
  filter(`0_B1` == 'unimodal symmetric') %>%
  mutate(is_robust = ifelse(`0_B1` == `1_B1`, T, F)) %>%
  dplyr::select(colnames(lhs_sets6), is_robust) %>% ungroup() %>%
  dplyr::select(-paramset) 

temp_class_for_tree$is_robust <- as.factor(temp_class_for_tree$is_robust)

temp.tree <- ctree(is_robust ~ basal_nitc_on_ratio + onbasalA1_off_ratio + A1_Aprime1_addon_ratio + A1_Aprime_prodon_ratio + r_addon_byA1_B1 , data = temp_class_for_tree)




