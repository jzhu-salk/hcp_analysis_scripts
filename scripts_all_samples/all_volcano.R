library(ggpubr)

# Figure 3C.
get_volcano <- function(deseq2_df, title){
  # plot genes in volcano plot
  volcano <- deseq2_df %>%
    filter(!is.na(padj)) %>%
    ggplot(mapping=aes(x=log2FoldChange, y=-log10(padj), color=Cutoff)) +
    geom_point(alpha=0.8, size=1) +
    scale_color_manual(values=c('dec' = 'blue', 'inc' = 'red', 'ns' = 'black', 'nfc' = 'black')) +
    scale_x_continuous(limits = c(0 - deseq2_df %>% pull(log2FoldChange) %>% abs %>% max,
                                  deseq2_df %>% pull(log2FoldChange) %>% abs %>% max)) +
    labs(title = title,
         subtitle = paste('Increasing Genes:', filter(deseq2_df, Cutoff == 'inc') %>% nrow,
                          '\nDecreasing Genes:', filter(deseq2_df, Cutoff == 'dec') %>% nrow,
                          '\nNS/NFC Genes:', filter(deseq2_df, Cutoff %in% c('nfc', 'ns')) %>% nrow),
         x=expression('Log'[2]*'(Fold Change)'),
         y=expression(-'Log'[10]*'(P-Adjusted Value)')) +
    geom_vline(xintercept=c(1, -1), linetype = 'dashed') +
    geom_hline(yintercept=-log10(0.05), linetype = 'dashed') +
    theme_minimal() +
    theme(legend.position = 'none',
          plot.subtitle = element_text(size = 9, hjust = 0),
          #panel.grid.minor.y = element_blank()
    )
  
  return(volcano)
}

# Individual v Rest ----
hcp_volcano = get_volcano(hcp_df, 'hCP vs Others Volcano')
hnp_volcano = get_volcano(hnp_df, 'hNP vs Others Volcano')
hnat_volcano = get_volcano(hnat_df, 'hNAT vs Others Volcano')

ggarrange(hcp_volcano, hnp_volcano, hnat_volcano, ncol = 3)

# Generate plot highlighting significant genes with most extreme fold change
extreme_fc = hcp_df %>% 
  filter(padj < 0.05) %>%
  arrange(log2FoldChange) %>%
  filter(row.names(.) %in% row.names(union(head(., 10), tail(., 10)))) %>%
  mutate(Symbol = coalesce(Symbol, ENSG))

highlight_volcano = 
  ggplot(mapping=aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data = hcp_df %>% filter(!is.na(padj)), aes(color = Cutoff), alpha=0.2, size = 1) + 
  geom_point(data = extreme_fc, alpha = .5, size = 2, color = 'green') +
  scale_color_manual(values=c('dec' = 'blue', 'inc' = 'red', 'ns' = 'black', 'nfc' = 'black')) +
  scale_x_continuous(limits = c(0 - hcp_df %>% pull(log2FoldChange) %>% abs %>% max,  
                                hcp_df %>% pull(log2FoldChange) %>% abs %>% max)) +
  labs(title = 'hCP vs Others Highlighted Genes',
       subtitle = paste('Increasing Genes:', filter(hcp_df, Cutoff == 'inc') %>% nrow,
                        '\nDecreasing Genes:', filter(hcp_df, Cutoff == 'dec') %>% nrow, 
                        '\nNS/NFC Genes:', filter(hcp_df, Cutoff %in% c('nfc', 'ns')) %>% nrow),
       x=expression('Log'[2]*'(Fold Change)'),
       y=expression(-'Log'[10]*'(P-Adjusted Value)')) +
  geom_vline(xintercept=c(1, -1), linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), linetype = 'dashed') +
  geom_label_repel(data = extreme_fc, aes(label=Symbol), direction = 'both', min.segment.length = 0) +
  theme_minimal() +
  theme(legend.position = 'none',
        plot.subtitle = element_text(size = 9, hjust = 0), 
        #panel.grid.minor.y = element_blank()
  )
plot(highlight_volcano)    

