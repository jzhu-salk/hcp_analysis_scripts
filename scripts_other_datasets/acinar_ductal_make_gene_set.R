acinar_ductal_gene_lists = list.files('gsea/acinar_ductal/', pattern = 'ranked_gene_list', recursive = T, full.names = T)

for (single_acinar_ductal_gene_list in acinar_ductal_gene_lists){
  print(single_acinar_ductal_gene_list)
  
  gene_list = read.csv(single_acinar_ductal_gene_list, sep='\t')
  
  gene_set_name = single_acinar_ductal_gene_list %>% basename %>% sub('ranked_gene_list_(.*)_\\d*[.]tsv', '\\1', .)
  print(gene_set_name)
  top_500 = gene_list %>% head(500) %>% pull(NAME)
  bot_500 = gene_list %>% tail(500) %>% pull(NAME)

  top_500_formatted = rep(paste0(gene_set_name, '_top_500'), 2) %>%
    c(., top_500)
  # write.table(paste(top_500_formatted, collapse = '\t'), file = paste0('gsea/acinar_ductal/', gene_set_name, '_top_500_gene_set.tsv'), sep = '\t', quote = F, row.names = F, col.names = F)

  bot_500_formatted = rep(paste0(gene_set_name, '_bot_500'), 2) %>%
    c(., bot_500)
  # write.table(paste(bot_500_formatted, collapse = '\t'), file = paste0('gsea/acinar_ductal/', gene_set_name, '_bot_500_gene_set.tsv'), sep = '\t', quote = F, row.names = F, col.names = F)
}

# Concatenate all output files in command line:
# for i in *tsv ; do cat $i >> acinar_ductal_gene_set.gmt ; done