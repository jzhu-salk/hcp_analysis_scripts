library(DESeq2)
library(tidyverse)
library(data.table)

# Conflict with "select" function, so detach in order to use "Reduce" function
detach("package:org.Hs.eg.db")
detach("package:AnnotationDbi")

# individually list each sample quantification file, read each file, then merge on gene_id
ht_gene_mat = list.files(paste0(ht_dir, '/align_quant'), pattern = '.*ReadsPerGene.out.tab', full.names=T) %>%
  # unstranded library prep, use col 2
  lapply(function(x) fread(file=x, skip = 4, col.names = c('gene_id',  sub('_ReadsPerGene.out.tab', '', basename(eval(x))), 'fwd', 'rev')) %>% select(c(1,2))) %>%
  # merge on gene_id column
  Reduce(function(x, y) merge(x, y, by='gene_id'), .) %>%
  column_to_rownames('gene_id')

# add 't_' because some samples have incompatible number headers
colnames(ht_gene_mat) = paste0('t_', colnames(ht_gene_mat))

# Metadata from GDC download

ht_metadata = read.csv(paste0(ht_dir, '/metadata/gdc_sample_sheet.2024-11-26.tsv'), sep='\t') %>% 
  filter(grepl('gene_counts', File.Name)) %>%
  select(where(~n_distinct(.) > 1)) %>%
  merge(read.csv(paste0(ht_dir, '/metadata/sample.tsv'), sep='\t') %>% select(where(~n_distinct(.) > 1)), by.x = c('Sample.ID', 'Case.ID'), by.y = c('sample_submitter_id', 'case_submitter_id')) %>%
  merge(read.csv(paste0(ht_dir, '/metadata/clinical.tsv'), sep='\t') %>% select(where(~n_distinct(.) > 1)), by.x = c('Case.ID', 'case_id'), by.y = c('case_submitter_id', 'case_id')) %>%
  merge(read.csv(paste0(ht_dir, '/metadata/phs001611.v1.pht009161.v1.p1.c1.Organoid_Profiling_PC_Sample_Attributes.GRU.txt'), sep='\t', skip=10) %>% select(where(~n_distinct(.) > 1)), by.x = 'Sample.ID', by.y = 'SAMPLE_ID')

ht_colData = ht_metadata %>%
  mutate(File.Name = paste0('t_', File.Name) %>% gsub('-.*', '', .),
         batch = 'hT',
         BODY_SITE_mod = sub(' .*', '', BODY_SITE),
         type = sub('(h.).*', '\\1', CORRESPONDING_ORGANOID)) %>%
  arrange(File.Name, colnames(ht_gene_mat)) %>% 
  filter(BODY_SITE_mod == 'Pancreas') %>%
  column_to_rownames('File.Name')

ht_gene_mat = ht_gene_mat %>% dplyr::select(row.names(ht_colData))

ht_dds = DESeqDataSetFromMatrix(countData = ht_gene_mat,
                                    colData = ht_colData,
                                    design = ~ tissue_type)
ht_dds %>% counts %>% rowSums %>% summary
# remove genes with low reads
has_reads = rowMeans(counts(ht_dds)) > 1
ht_dds_filtered = ht_dds[has_reads,]
ht_dds_filtered %>% counts %>% rowSums %>% summary
ht_dds_final = DESeq(ht_dds_filtered)

