library(tidyverse)
library(ggrepel)
library(magrittr)
library(biomaRt)
# library(GenomicFeatures)

setwd('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/GSE92872/')

plotdir <- '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/'
# load and reformat bulk data
bulkcounts <- as_tibble(read.csv(file = 'GSE92872_CROP-seq_Jurkat_TCR.count_matrix.csv', stringsAsFactors = F))

bulk_coldata <- as.data.frame(bulkcounts[1:3,])
rownames(bulk_coldata) <- bulk_coldata$sample_name
bulk_coldata_tall <- t(bulk_coldata)[2:ncol(bulk_coldata),]

bulk_rowdata <- as.data.frame(bulkcounts[5:nrow(bulkcounts),1])
rownames(bulk_rowdata) <- bulk_rowdata$sample_name
bulk_rowdata %<>% dplyr::rename(gene_name = sample_name)

bulk_countmat <- as.data.frame(bulkcounts[5:nrow(bulkcounts),])
rownames(bulk_countmat) <- bulk_countmat$sample_name
bulk_countmat %<>% dplyr::select(-sample_name) %>% as.matrix()

target_genes <- unique(as.character(as.matrix(bulk_coldata[2,2:ncol(bulk_coldata)]))) #includes CTRL (control)


# load CROP-seq data
# sccounts <- as_tibble(read.csv(file = 'GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv', stringsAsFactors = F))

# pull list of paralogs from ensembl
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # version 105, Dec 2021
human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105) # version 105, Dec 2021
geneParaList <- getBM(attributes = c("ensembl_gene_id", 
                                     "external_gene_name",
                                     "hsapiens_paralog_ensembl_gene", 
                                     "hsapiens_paralog_associated_gene_name",
                                     'hsapiens_paralog_subtype',
                                     'hsapiens_paralog_orthology_type',
                                     'hsapiens_paralog_perc_id',
                                     'hsapiens_paralog_perc_id_r1'),
                      filters = 'external_gene_name',
                      values = target_genes,
                      mart = human)

# load "gene lengths", (i.e., union of all Refseq transcript exons/UTRs per gene) from hg38 genes
lengthtbl<- as_tibble(read.csv(file = '~/code/grn_nitc/Resources/hg38.ncbiRefSeq.txLengthPerGene.csv', header = T, stringsAsFactors = F)) %>%
  dplyr::rename(gene_name = gene_id)

# look for nitc in bulk results
bulk_counts_tall <- as.data.frame(bulk_countmat) %>%
  mutate(gene_name = rownames(bulk_countmat)) %>%
  as_tibble() %>%
  pivot_longer(cols = -gene_name, names_to = 'sampleID', values_to = 'count') 

bulk_counts_tall$count <- as.numeric(as.character(bulk_counts_tall$count))

bulk_totalcounts <- bulk_counts_tall %>%
  group_by(sampleID) %>%
  summarise(total_counts = sum(count))

countplot <- bulk_totalcounts %>% dplyr::select(total_counts) %>% unique() %>% ggplot() + geom_histogram(aes(total_counts)) + theme_classic() + ggtitle('total counts per bulk sample')
ggsave(countplot, file = paste0(plotdir, 'total_counts_per_bulk_sample.pdf'))

bulk_counts_tall %<>% 
  inner_join(bulk_totalcounts, by = 'sampleID') %>%
  mutate(RPM = count*1000000/total_counts) %>%
  inner_join(as.data.frame(bulk_coldata_tall) %>% 
               mutate(sampleID = rownames(bulk_coldata_tall)), by = 'sampleID') %>%
  filter(total_counts > 1e6)

# absolute change in TPM after KO per paralog
bulk_counts_tall %>%
  inner_join(lengthtbl, by = 'gene_name') %>%
  mutate(RPK = count*1e3/length) %>%
  group_by(gene) %>%
  mutate(sumRPKpm = sum(RPK)/1e6,
         TPM = RPK/sumRPKpm) %>%
  dplyr::select(-c('RPK','sumRPKpm','length'))
         
pseud = 1
bulk_FC_perTarget <- list()
bulk_FC_perTarget_perParalog <- list()
for(crispr_target in target_genes) {
  
  paralogs <- geneParaList %>%
    filter(external_gene_name == crispr_target)
  
  bulk_sub <- bulk_counts_tall %>%
    filter(gene %in% c(crispr_target, 'CTRL'),
           gene_name %in% c(crispr_target, paralogs$hsapiens_paralog_associated_gene_name))
  
  bulk_sub_sum_rpm <- bulk_sub %>%
    group_by(gene_name, gene, condition) %>%
    summarise(meanRPM = mean(RPM + pseud),
              sdRPM = sd(RPM + pseud),
              semRPM = sdRPM/sqrt(length(RPM)),
              CRISPR_target = crispr_target) %>%
    group_by(CRISPR_target, gene_name, condition) %>%
    pivot_wider(names_from = gene, values_from = c(meanRPM, sdRPM, semRPM)) %>%
    mutate(meanFC = eval(as.symbol(paste0('meanRPM_', crispr_target)))/meanRPM_CTRL,
           sdFC = meanFC*sqrt((eval(as.symbol(paste0('sdRPM_', crispr_target)))/eval(as.symbol(paste0('meanRPM_', crispr_target))))^2 + 
                                (sdRPM_CTRL/meanRPM_CTRL)^2),
           lfc_RPM = log2(meanFC),
           lfc_RPM_up = log2(meanFC+sdFC),
           lfc_RPM_dn = log2(max((meanFC-sdFC),0.001)))
  
  bulk_summary_FC_rpm <- bulk_sub_sum_rpm %>%
    group_by(condition) %>%
    filter(gene_name != crispr_target) %>%
    summarise(gmeanParalogFC = exp(mean(log(meanFC))),
              log2_gmeanParalogFC = log2(gmeanParalogFC),
              nParalogs = length(meanFC)) %>%
    mutate(CRISPR_target = crispr_target)
  
  p1 <- ggplot(bulk_sub, aes(gene, RPM)) +
    facet_grid(condition~gene_name) +
    geom_point(aes(color = gene)) +
    theme_bw() +
    ggtitle(paste0('bulk gene expression after ', crispr_target, ' knockout'))
  
  ggsave(p1, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO.pdf'), height = 4, width = 2+nrow(paralogs))
  p2 <- ggplot(bulk_sub, aes(gene, log2(RPM+1))) +
    facet_grid(condition~gene_name) +
    geom_point(aes(color = gene)) +
    theme_bw() +
    ylab('log2(RPM+1)') +
    ggtitle(paste0('bulk gene expression after ', crispr_target, ' knockout'))
  
  ggsave(p2, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log.pdf'), height = 4, width = 2+nrow(paralogs))
  
  p3 <- ggplot() +
    geom_point(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), lfc_RPM)) +
    geom_errorbar(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), ymin = lfc_RPM_dn, ymax = lfc_RPM_up)) +
    theme_classic() +
    geom_text_repel(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), lfc_RPM, label = gene_name)) +
    facet_wrap(~condition) +
    ggtitle(paste0('bulk gene expression change relative to control\nafter ', crispr_target, ' knockout')) +
    ylab('log2(fold change) relative to control') +
    xlab('log2(meanRPM+1) in controls')
  
  ggsave(p3, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_FC.pdf'), height = 4, width = 8)
  
  p4 <- ggplot() +
    geom_point(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), lfc_RPM)) +
    # geom_errorbar(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), ymin = lfc_RPM_dn, ymax = lfc_RPM_up)) +
    theme_classic() +
    geom_text_repel(data = bulk_sub_sum, aes(log2(meanRPM_CTRL+1), lfc_RPM, label = gene_name)) +
    facet_wrap(~condition) +
    ggtitle(paste0('bulk gene expression change relative to control\nafter ', crispr_target, ' knockout')) +
    ylab('log2(fold change) relative to control') +
    xlab('log2(meanRPM+1) in controls')
  
  ggsave(p4, file = paste0('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_',crispr_target,'_KO_log_FCnoEB.pdf'), height = 4, width = 8)
  
  if(is.null(dim(bulk_FC_perTarget))) {
    bulk_FC_perTarget <- bulk_summary_FC
  } else {
    bulk_FC_perTarget %<>% bind_rows(bulk_summary_FC)
  }
  
  if(is.null(dim(bulk_FC_perTarget_perParalog))) {
    bulk_FC_perTarget_perParalog <- bulk_sub_sum
  } else {
    bulk_FC_perTarget_perParalog %<>% bind_rows(bulk_sub_sum)
  }
  
}

bulk_FC_EGR1_stim_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_text_repel(data = bulk_FC_perTarget_perParalog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, lfc_RPM, label = gene_name), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget %>% filter(CRISPR_target == 'EGR1', condition == 'stimulated'), aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_EGR1_stim_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_EGR1_stim_KO_FCplot.pdf', height = 4, width = 2)

bulk_FC_EGR1_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1'), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_text_repel(data = bulk_FC_perTarget_perParalog %>% filter(CRISPR_target != gene_name, CRISPR_target == 'EGR1'), aes(CRISPR_target, lfc_RPM, label = gene_name), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget %>% filter(CRISPR_target == 'EGR1'), aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_EGR1_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_EGR1_KO_FCplot.pdf', height = 7, width = 2)

bulk_FC_perTarget_FCplot <- ggplot() + 
  geom_point(data = bulk_FC_perTarget_perParalog %>% filter(CRISPR_target != gene_name), aes(CRISPR_target, lfc_RPM), position = position_jitter(seed = 8272, width = 0.1, height = 0)) +
  geom_bar(data = bulk_FC_perTarget, aes(CRISPR_target, log2_gmeanParalogFC), stat = 'identity', alpha = 0.3) +
  facet_grid(condition~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30)) +
  ylab('Log2(fold-change) after KO in points\nLog2(Geometric mean fold change) in bar') +
  xlab('CRISPR target') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nDoes not include reference gene change')
ggsave(bulk_FC_perTarget_FCplot, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulkExp_allParalogs_KO_log_FCnoEB.pdf', height = 7, width = 7)

paralogs %>% as_tibble()


bulk_FCvsPercID_EGR1_stim_perTarget <- ggplot(bulk_FC_perTarget_perParalog %>% 
                                                filter(CRISPR_target != gene_name,
                                                       CRISPR_target == 'EGR1',
                                                       condition == 'stimulated') %>%
                                                inner_join(geneParaList %>%
                                                             dplyr::rename(CRISPR_target = external_gene_name,
                                                                           gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
                                              aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  geom_point() +
  geom_text_repel(aes(label = gene_name)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene')
ggsave(bulk_FCvsPercID_EGR1_stim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_EGR1_stim_perTarget.pdf', height = 4, width = 4)

bulk_FCvsPercID_stim_perTarget <- ggplot() +
  geom_point(data = bulk_FC_perTarget_perParalog %>% 
               filter(CRISPR_target != gene_name,
                      condition == 'stimulated') %>%
               inner_join(geneParaList %>%
                            dplyr::rename(CRISPR_target = external_gene_name,
                                          gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
             aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene\nStimulated condition')
ggsave(bulk_FCvsPercID_stim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_stim_perTarget.pdf', height = 7, width = 7)

bulk_FCvsPercID_unstim_perTarget <- ggplot() +
  geom_point(data = bulk_FC_perTarget_perParalog %>% 
               filter(CRISPR_target != gene_name,
                      condition == 'unstimulated') %>%
               inner_join(geneParaList %>%
                            dplyr::rename(CRISPR_target = external_gene_name,
                                          gene_name = hsapiens_paralog_associated_gene_name), by = c('CRISPR_target', 'gene_name')), 
             aes(hsapiens_paralog_perc_id, lfc_RPM)) +
  facet_wrap(~CRISPR_target) +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab('Log2(fold-change) after KO') +
  ggtitle('Change in paralog expression (RPM+1) after reference gene KO\nVs. Percent sequence similarity to reference gene\nUnstimulated condition')
ggsave(bulk_FCvsPercID_unstim_perTarget, file = '/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/exploration/bulk_FCvsPercID_unstim_perTarget.pdf', height = 7, width = 7)

