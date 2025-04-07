# # All Samples ----
# DDS
source('scripts_all_samples/all_dds.R')
# PCA Plot
source('scripts_all_samples/all_pca.R')
# Heat Map
source('scripts_all_samples/all_heatmap.R')
# DEG
source('scripts_all_samples/all_deg.R')
# Volcano Plot
source('scripts_all_samples/all_volcano.R')

# # hCP Samples ----
# DDS
source('scripts_hcp_samples/hcp_dds.R')
# PCA Plot
source('scripts_hcp_samples/hcp_pca.R')
# NMF
source('scripts_hcp_samples/hcp_nmf.R')
# NMF-Based GSEA
source('scripts_hcp_samples/hcp_gsea.R')

# # Published Dataset Samples ----
# Acinar Ductal Samples DDS
source('scripts_other_datasets/acinar_ductal_dds.R')
# Acinar Ductal Samples Create Gene Sets
source('scripts_other_datasets/acinar_ductal_make_gene_set.R')
# Acinar Ductal Samples Run ssGSEA
source('scripts_other_datasets/acinar_ductal_ssgsea.R')

# hT DDS
source('scripts_other_datasets/ht_dds.R')
# hT DEG
source('scripts_other_datasets/ht_deg.R')
# hT vs hCP Comparison
source('scripts_other_datasets/ht_hcp_compare.R')

# save.image(file = 'hcp_github.RData', compress = 'gzip')
# > sessionInfo()
# R version 4.4.3 (2025-02-28)
# Platform: x86_64-redhat-linux-gnu
# Running under: Fedora Linux 40 (Server Edition)
# 
# Matrix products: default
# BLAS/LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: America/Los_Angeles
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] org.Hs.eg.db_3.20.0         AnnotationDbi_1.68.0        doParallel_1.0.17           iterators_1.0.14            foreach_1.5.2              
# [6] NMF_0.25                    cluster_2.1.8               rngtools_1.5.2              registry_0.5-1              ggpubr_0.6.0               
# [11] pheatmap_1.0.12             colorBlindness_0.1.9        ggrepel_0.9.6               openxlsx_4.2.7.1            data.table_1.16.4          
# [16] lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                
# [21] readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0            
# [26] DESeq2_1.46.0               SummarizedExperiment_1.36.0 Biobase_2.66.0              MatrixGenerics_1.18.1       matrixStats_1.5.0          
# [31] GenomicRanges_1.58.0        GenomeInfoDb_1.42.1         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0        
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3               rlang_1.1.5             magrittr_2.0.3          gridBase_0.4-7          compiler_4.4.3         
# [6] RSQLite_2.3.9           reshape2_1.4.4          png_0.1-8               vctrs_0.6.5             pkgconfig_2.0.3        
# [11] crayon_1.5.3            fastmap_1.2.0           backports_1.5.0         XVector_0.46.0          labeling_0.4.3         
# [16] tzdb_0.4.0              UCSC.utils_1.2.0        bit_4.5.0.1             zlibbioc_1.52.0         cachem_1.1.0           
# [21] jsonlite_1.8.9          blob_1.2.4              DelayedArray_0.32.0     BiocParallel_1.40.0     broom_1.0.7            
# [26] R6_2.5.1                stringi_1.8.4           RColorBrewer_1.1-3      car_3.1-3               Rcpp_1.0.14            
# [31] Matrix_1.7-2            timechange_0.3.0        tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
# [36] codetools_0.2-20        plyr_1.8.9              lattice_0.22-6          withr_3.0.2             KEGGREST_1.46.0        
# [41] gridGraphics_0.5-1      zip_2.3.1               Biostrings_2.74.1       BiocManager_1.30.25     pillar_1.10.1          
# [46] carData_3.0-5           generics_0.1.3          hms_1.1.3               munsell_0.5.1           scales_1.3.0           
# [51] glue_1.8.0              tools_4.4.3             locfit_1.5-9.10         ggsignif_0.6.4          cowplot_1.1.3          
# [56] grid_4.4.3              colorspace_2.1-1        GenomeInfoDbData_1.2.13 Formula_1.2-5           cli_3.6.3              
# [61] S4Arrays_1.6.0          gtable_0.3.6            rstatix_0.7.2           digest_0.6.37           SparseArray_1.6.1      
# [66] farver_2.1.2            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7              bit64_4.6.0-1  