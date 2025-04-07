hcp_gsea_files = list.dirs('gsea') %>%
  grep('vov_0693', ., value = T) %>%
  list.files(., pattern = '*report.*tsv', recursive = T, full.names = T)

hcp_gsea_df = lapply(hcp_gsea_files, function(y){
    read.csv(y, sep='\t') %>% 
      dplyr::select(c(NAME, SIZE, NES, FDR.q.val)) %>%
      mutate(Basis = gsub('.*0693_(\\d)_vs.*', '\\1', y)[[1]])
  }) %>%
  bind_rows

filtered_hcp_gsea_df = hcp_gsea_df %>%
  group_by(NAME) %>%
  filter(any(FDR.q.val < .05)) %>%
  ungroup %>%
  mutate(NAME = sub('HALLMARK_', '', NAME),
         Significance = ifelse(FDR.q.val < .05, 'Sig', 'NS'))

# Figure S4B.
hcp_basis_gsea = filtered_hcp_gsea_df %>%
  arrange(NAME) %>%
  ggplot() +
  geom_point(aes(x = NES, y = NAME, size = SIZE, color = Basis, shape = Significance)) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual(values = c('1' = "#0072B2", '2' = "#D55E00", '3' = "#CC79A7")) +
  theme_minimal() +
  labs(title = 'GSEA hCP Basis Comparison',
       y = element_blank(),
       size = 'Gene Set Size',
       color = 'Basis') +
  theme(strip.background = element_rect(colour = "black"),
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text.y = element_text(face = 'bold'))

hcp_basis_gsea
