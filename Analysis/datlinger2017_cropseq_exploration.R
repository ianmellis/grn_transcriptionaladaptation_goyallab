library(tidyverse)
library(magrittr)
library(biomaRt)

setwd('/Volumes/IAMYG1/grn_nitc_data/CROP-seq/Datlinger2017/GSE92872/')

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
sccounts <- as_tibble(read.csv(file = 'GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv', stringsAsFactors = F))

# pull list of paralogs from ensembl
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # version 105, Dec 2021
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

# look for nitc in bulk results
bulk_counts_tall <- as.data.frame(bulk_countmat) %>%
  mutate(gene_name = rownames(bulk_countmat)) %>%
  as_tibble() %>%
  pivot_longer(cols = -gene_name, names_to = 'sampleID', values_to = 'count') 

bulk_counts_tall$count <- as.numeric(as.character(bulk_counts_tall$count))

bulk_totalcounts <- bulk_counts_tall %>%
  group_by(sampleID) %>%
  summarise(total_counts = sum(count))

bulk_counts_tall %<>% 
  inner_join(bulk_totalcounts, by = 'sampleID') %>%
  mutate(RPM = count*1000000/total_counts) %>%
  inner_join(as.data.frame(bulk_coldata_tall) %>% mutate(sampleID = rownames(bulk_coldata_tall)), by = 'sampleID')

for(crispr_target in target_genes) {
  
  paralogs <- geneParaList %>%
    filter(external_gene_name == crispr_target)
  
  bulk_sub <- bulk_counts_tall %>%
    filter(gene %in% c(crispr_target, 'CTRL'),
           gene_name %in% c(crispr_target, paralogs$hsapiens_paralog_associated_gene_name))
  
  bulk_sub_sum <- bulk_sub %>%
    group_by(gene_name, gene, condition) %>%
    summarise(meanRPM = mean(RPM),
              sdRPM = sd(RPM),
              semRPM = sd(RPM)/sqrt(length(RPM))) %>%
    group_by(gene_name, condition) %>%
    pivot_wider(names_from = gene, values_from = c(meanRPM, sdRPM, semRPM)) %>%
    mutate(meanFC = eval(as.symbol(paste0('meanRPM_', crispr_target)))/meanRPM_CTRL,
           sdFC = meanFC*sqrt((eval(as.symbol(paste0('sdRPM_', crispr_target)))/eval(as.symbol(paste0('meanRPM_', crispr_target))))^2 + 
                                (sdRPM_CTRL/meanRPM_CTRL)^2),
           lfc_RPM = log2(meanFC),
           lfc_RPM_up = log2(meanFC+sdFC),
           lfc_RPM_dn = log2(max((meanFC-sdFC),0.001)))
  
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
  
}
