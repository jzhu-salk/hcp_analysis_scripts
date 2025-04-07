library(DESeq2)
library(tidyverse)
library(data.table)
library(openxlsx)

# Read from STAR quantification output, column 4 (reverse stranded library prep)
gene_mat = list.files('../align_quant/', pattern = 'ReadsPerGene.out.tab', full.names = T) %>%
  lapply(function(x) fread(file=x, sep='\t', skip=4, col.names=c('genes', 'unstranded', 'forward', gsub('_ReadsPerGene.out.tab', '', basename(eval(x))))) %>% select(c(1,4))) %>%
  Reduce(function(x, y) merge(x, y, by='genes'), .) %>%
  column_to_rownames('genes')

# Error in naming, rename hNAT files
gene_mat = gene_mat %>%
  rename_with(.fn = ~sub('02', '01', .), .cols = grep('^hNAT02', colnames(gene_mat))) %>%
  rename_with(.fn = ~sub('04', '03', .), .cols = grep('^hNAT04', colnames(gene_mat))) %>%
  rename_with(.fn = ~sub('07', '06', .), .cols = grep('^hNAT07', colnames(gene_mat))) %>%
  rename_with(.fn = ~sub('09', '08', .), .cols = grep('^hNAT09', colnames(gene_mat)))

colnames(gene_mat) = colnames(gene_mat) %>%
  sub('0(\\d)', '\\1', .) %>%
  sub("_(.).*", "\\1", .)

colData = colnames(gene_mat) %>%
  data.frame(row.names = .,
             condition = sub('\\d.*', '', .))

dds = DESeqDataSetFromMatrix(countData = gene_mat,
                             colData = colData,
                             design = ~ condition)
dds %>% counts %>% rowSums %>% summary
min_group_size = colData(dds)['condition'] %>% table %>% min
# remove genes with low reads
has_reads = rowSums(counts(dds) > 10 ) > min_group_size
dds_filtered = dds[has_reads,]
dds_filtered %>% counts %>% rowSums %>% summary
dds_final = DESeq(dds_filtered)

# Prepare Normalized Counts
library(org.Hs.eg.db)
all_samples_norm_counts = dds_final %>% counts(normalized=T) %>% data.frame %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, keys = gsub('\\..*', '', rownames(.)), column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
# write.table(x = all_samples_norm_counts, file = 'norm_counts/all_samples_norm_counts.tsv', quote = F, sep = '\t')

gsea_counts = all_samples_norm_counts %>%
  filter(!is.na(Symbol)) %>%
  remove_rownames %>%
  distinct(Symbol, .keep_all = T) %>%
  column_to_rownames('Symbol') %>%
  mutate(DESCRIPTION = NA) %>%
  relocate(DESCRIPTION) %>%
  rownames_to_column(var = 'NAME')
# write.table(x = gsea_counts, file = 'gsea/all_samples_gsea_norm_counts.tsv', quote = F, row.names = F, sep='\t')
# Will need to manually add additional lines to top of this file and rename as gct for GSEA:
# #1.2
# 14560 37