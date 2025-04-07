hcp_df # hCP vs hNAT/hNP
ht_df # Tumor vs Normal

merged_df = merge(ht_df  %>% dplyr::select(ENSG, Symbol, Cutoff, log2FoldChange),
                  hcp_df %>% dplyr::select(ENSG, Symbol, Cutoff, log2FoldChange),
                  by = c('ENSG', 'Symbol'), suffixes = c('_ht', '_hcp')) %>%
  filter(Cutoff_ht != 'ns' & Cutoff_hcp != 'ns') %>%
  mutate(merged_cutoff = paste0('hT_', Cutoff_ht, '_hCP_', Cutoff_hcp),
         Symbol = coalesce(Symbol, ENSG),
         ht_hcp_lfc_magnitude_diff = log2FoldChange_ht - log2FoldChange_hcp)

top_magnitude_diff = merged_df %>% filter(ENSG %in% 
                                            union(merged_df %>% arrange(ht_hcp_lfc_magnitude_diff) %>% head(4) %>% pull(ENSG),
                                                  merged_df %>% arrange(ht_hcp_lfc_magnitude_diff) %>% tail(4) %>% pull(ENSG)))

# Figure S4E.
hcp_ht_comparison_plot = ggplot(mapping = aes(x = log2FoldChange_ht, y = log2FoldChange_hcp)) +
  geom_point(data = merged_df) +
  geom_label_repel(data = top_magnitude_diff, mapping = aes(label = Symbol), min.segment.length = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
  geom_hline(yintercept = c(-1, 1), linetype = 'dashed') +
  labs(title = 'Fold Change Comparison of Significant Genes',
       x = expression('Tumor vs Normal Log'[2]*'(Fold Change)'),
       y = expression('hCP vs hNAT/hNP Log'[2]*'(Fold Change)')) +
  theme_minimal()
plot(hcp_ht_comparison_plot)
