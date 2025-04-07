# Read in STAR quant outputs
acinar_ductal_gene_mat = list.files(paste0(acinar_ductal_dir, '/align_quant'), pattern = 'ReadsPerGene.out.tab', full.names = T) %>%
  lapply(function(x) fread(file=x, sep='\t', skip=4, col.names=c('genes', 'unstranded', 'forward', gsub('_ReadsPerGene.out.tab', '', basename(eval(x))))) %>% dplyr::select(c(1,4))) %>%
  Reduce(function(x, y) merge(x, y, by='genes'), .) %>%
  column_to_rownames('genes')

acinar_ductal_colData = read.csv(paste0(acinar_ductal_dir, 'SraRunTable.csv')) %>%
  dplyr::select(where(~n_distinct(.) > 1)) %>%
  column_to_rownames('Run') %>%
  mutate(condition = sub(' ', '_', condition),
         cell_type_genotype_condition = paste0(cell_type, '_', genotype, '_', condition),
         genotype_condition = paste0(genotype, '_', condition))

acinar_ductal_dds = DESeqDataSetFromMatrix(countData = acinar_ductal_gene_mat,
                                           colData = acinar_ductal_colData,
                                           design = ~ 1) 

acinar_ductal_dds %>% counts %>% rowSums %>% summary
# remove genes with low reads
acinar_ductal_has_reads =  rowMeans(acinar_ductal_dds %>% counts) > 1

acinar_ductal_dds_filtered = acinar_ductal_dds[acinar_ductal_has_reads,]
acinar_ductal_dds_filtered %>% counts %>% rowSums %>% summary
acinar_ductal_dds_final = DESeq(acinar_ductal_dds_filtered)


# Writing Out Files ----
# Normalized Counts
library(org.Hs.eg.db)
acinar_ductal_norm_counts = acinar_ductal_dds_final %>% counts(normalized=T) %>% data.frame %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, keys = gsub('\\..*', '', rownames(.)), column="SYMBOL", keytype="ENSEMBL", multiVals="first"))

acinar_ductal_gsea_counts = acinar_ductal_norm_counts %>%
  filter(!is.na(Symbol)) %>%
  remove_rownames %>%
  distinct(Symbol, .keep_all = T) %>%
  column_to_rownames('Symbol') %>%
  mutate(DESCRIPTION = NA) %>%
  relocate(DESCRIPTION) %>%
  rownames_to_column(var = 'NAME')

# write.table(acinar_ductal_gsea_counts, file = 'gsea/acinar_ductal_gsea_counts.txt', sep = '\t', row.names = F, quote = F)
