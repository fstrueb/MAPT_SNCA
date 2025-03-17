MSBB data anlysis pertaining to Figure 4
================
<felix.struebing@med.uni-muenchen.de>
2025-03-17

- [0.1 Proteome](#01-proteome)
- [0.2 RNA](#02-rna)
- [0.3 Session info](#03-session-info)

``` r
library(tidyverse)
library(ggraph)
library(patchwork)
library(RColorBrewer)
library(viridis)
set.seed(42)
```

## 0.1 Proteome

Get MSBB proteome data:

``` r
proteome <- read_tsv('/earth/public_data/datasets/MSBB/Proteomics/MSSM_Proteomics_PFC_PROTEINOUTPUT.txt')
metad <- read_csv('/earth/public_data/datasets/MSBB/Proteomics/MSBB_Proteomics_PFC_TRAITS.csv') |> 
    separate(RunName, into = c('blabla', 'ID'), sep = '_RAW_') |> 
    separate(ID, into = c('bx', 'dig4', 'dig2'), sep = '_') |> 
    mutate(name = paste(bx, dig2, dig4, sep = '_'))
fasta_names <- Biostrings::readDNAStringSet('/earth/public_data/datasets/MSBB/Proteomics/MSSM_Proteomics_PFC_DATABASE.fasta')
```

Distribution of BB stages:

``` r
metad |> 
    ungroup() |> 
    mutate(bbscore = as.factor(bbscore)) |> 
    dplyr::filter(!is.na(bbscore)) |> 
    ggplot(aes(x = bbscore, fill = bbscore)) +
    geom_bar(stat = 'count', color = 'black') +
    coord_flip() +
    scale_fill_brewer(guide = 'none', palette = 'YlOrRd')
```

![](gh_04_MSBB_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
proteome_anno <- proteome %>% 
    dplyr::filter(is.na(`Potential contaminant`)) %>% # remove contaminants
    dplyr::filter(is.na(Reverse)) %>% # remove peptides derived from the reversed part of the decoy database
    dplyr::select(2, `Mol. weight [kDa]`, starts_with('Intensity')) |> 
    dplyr::select(-Intensity) %>% 
    dplyr::mutate(uniprot = str_extract_all(`Majority protein IDs`, pattern = '(?<=\\|).*?(?=_HUMAN)'))
temp_l <- sapply(proteome_anno$uniprot, length)
which(temp_l == 0)
```

    ## [1] 65

``` r
proteome_anno[65,]$uniprot
```

    ## [[1]]
    ## character(0)

``` r
proteome_anno <- proteome_anno %>% 
    dplyr::filter(`Majority protein IDs` != 'sp|ApoE4|LGADMEDVR') %>% # single outlier that didn't return a regex match
    dplyr::mutate(uniprot_pasted = unlist(lapply(uniprot, function(x) {
        unlist(strsplit(x, split = '\\|')[1])[1]}))) %>% 
    dplyr::mutate(uniprot = str_remove(uniprot_pasted, pattern = '\\-[0-9]'))
# proteome_anno$uniprot
# saveRDS(proteome_anno$uniprot, 'export/uniprot_ids_list.rds')
temp_ids <- lapply(proteome_anno$uniprot, function(x) {
    unlist(strsplit(x, split = '\\|')[1])[1]
}) %>% unlist(.) %>% str_remove(., pattern = '\\-[0-9]')
length(temp_ids) == length(proteome_anno$uniprot)
```

    ## [1] TRUE

``` r
source('../map_uniprot_ids.R') # result is called resultTable
proteome_anno <- proteome_anno %>% 
    left_join(., resultsTable, by = c('uniprot' = 'From')) %>% 
    dplyr::select(uniprot_pasted, gene_symbol = To, everything(), -uniprot)
temp_ids <- which(grepl('Intensity', names(proteome_anno)))
names(proteome_anno)[temp_ids] <- str_sub(names(proteome_anno)[temp_ids], start = 11L)
# names(proteome_anno)
proteome_anno <- proteome_anno %>% 
    dplyr::mutate(gene_symbol = ifelse(is.na(gene_symbol), paste0('uniprot_', tolower(uniprot_pasted)), gene_symbol)) # if no gene symbol is identified, use lower case uniprot ID
```

Make gene symbol names unique by appending their UNIPROT id in case they
have the same basename:

``` r
table(duplicated(proteome_anno$gene_symbol))
```

    ## 
    ## FALSE  TRUE 
    ##  5301   749

``` r
proteome_anno <- proteome_anno %>% 
    dplyr::mutate(gene_symbol_unique = ifelse(duplicated(gene_symbol), paste0(uniprot_pasted, '_', gene_symbol), gene_symbol))
table(duplicated(proteome_anno$gene_symbol_unique))
```

    ## 
    ## FALSE 
    ##  6050

``` r
proteome_mat <- proteome_anno %>% 
    dplyr::select(-uniprot_pasted, -`Mol. weight [kDa]`, -`Majority protein IDs`, -gene_symbol) %>% 
    column_to_rownames('gene_symbol_unique') %>% 
    as.matrix(.)
# names(proteome_mat)
# head(proteome_mat)
```

``` r
library(ggpubr)
proteome_cor <- proteome_anno %>% 
    dplyr::filter(grepl('SNCA|MAPT', gene_symbol)) %>% 
    dplyr::select(gene_symbol_unique, starts_with('b')) %>% 
    pivot_longer(-gene_symbol_unique) %>% 
    pivot_wider(names_from = 'gene_symbol_unique', values_from = 'value') %>% 
    left_join(metad, by = 'name') %>% 
    mutate(bbscore = as.factor(bbscore)) |> 
    dplyr::filter(!is.na(bbscore))
```

``` r
a <- proteome_cor |> 
    ggplot(aes(x = log10(`P10636-2_MAPT`), y = log10(SNCA))) +
    geom_point(alpha = 0.8, aes(color = bbscore)) +
    geom_smooth(method = 'lm', color = 'dodgerblue1') +
    stat_cor(label.y.npc = 'bottom') +
    labs(y = 'Alpha-Synuclein', x = '0N3R') +
    scale_color_brewer(guide = 'none', palette = 'YlOrRd') +
    theme(axis.text.x = element_text(angle = 90))
b <- proteome_cor |> 
    ggplot(aes(x = log10(`P10636-6_MAPT`), y = log10(SNCA))) +
    geom_point(alpha = 0.8, aes(color = bbscore)) +
    geom_smooth(method = 'lm', color = 'dodgerblue1') +
    stat_cor(label.y.npc = 'bottom') +
    labs(x = '0N4R', y = '') +
    scale_color_brewer(guide = 'none', palette = 'YlOrRd') +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90))
c <- proteome_cor |> 
    ggplot(aes(x = log10(`P10636-7_MAPT`), y = log10(SNCA))) +
    geom_point(alpha = 0.8, aes(color = bbscore)) +
    geom_smooth(method = 'lm', color = 'dodgerblue1') +
    stat_cor(label.y.npc = 'bottom') +
    labs(x = '1N4R', y = '') +
    scale_color_brewer(guide = 'none', palette = 'YlOrRd') +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90))
d <- proteome_cor |> 
    ggplot(aes(x = log10(`P10636-8_MAPT`), y = log10(SNCA))) +
    geom_point(alpha = 0.8, aes(color = bbscore)) +
    geom_smooth(method = 'lm', color = 'dodgerblue1') +
    stat_cor(label.y.npc = 'bottom') +
    labs(x = '2N4R', y = '') +
    scale_color_brewer(guide = 'none', palette = 'YlOrRd') +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90))
patchwork::wrap_plots(a, b, c, d, nrow = 1)
```

![](gh_04_MSBB_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## 0.2 RNA

Get RNA-seq identifiers:

``` r
metad_rna <- read_tsv('/earth/public_data/datasets/MSBB/RNAseq/bamfiles/SYNAPSE_METADATA_MANIFEST.tsv') |> 
    filter(individualID %in% metad$individualIdentifier) |> 
    filter(grepl('sort.coord.bam', name))
```

``` r
table(duplicated(metad_rna$individualID))
```

    ## 
    ## FALSE  TRUE 
    ##   223   670

``` r
metad_rna |> dplyr::count(tissue)
```

    ## # A tibble: 5 × 2
    ##   tissue                      n
    ##   <chr>                   <int>
    ## 1 frontal pole              236
    ## 2 inferior frontal gyrus    205
    ## 3 parahippocampal gyrus     191
    ## 4 superior temporal gyrus   214
    ## 5 <NA>                       47

``` r
metad_rna_fil <- metad_rna |> 
    filter(tissue == 'superior temporal gyrus') |>
    # filter(!grepl('hB', name)) |> 
    add_count(individualID) |> arrange(-n)
rna_files <- data.frame(bampath = list.files('/earth/public_data/datasets/MSBB/RNAseq/bamfiles', pattern = 'sort.coord.bam$', full.names = TRUE)) |> 
    mutate(basename = basename(bampath))
rna_files$size = sapply(rna_files$bampath, file.size, USE.NAMES = FALSE)
metad_rna_fil <- left_join(metad_rna_fil, rna_files, by = c('name' = 'basename')) |> 
    group_by(individualID) |> 
    # slice_max(size, n = 1) %>%
    dplyr::filter(!is.na(bampath))
table(duplicated(metad_rna_fil$individualID))
```

    ## 
    ## FALSE  TRUE 
    ##   100     1

``` r
table(is.na(metad_rna_fil$bampath))
```

    ## 
    ## FALSE 
    ##   101

Use featureCounts for RNA feature counting:

``` r
library(Rsubread)
# mycounts <- featureCounts(files = metad_rna_fil$bampath,
#                           annot.inbuilt = 'hg19',
#                           isPairedEnd = FALSE,
#                           nthreads = 18
# )
# saveRDS(mycounts, 'export/msbb_rnaseq_proteome_matchedCohort_featureCounts_supTempGyrus.rds')
mycounts <- readRDS('../export/msbb_rnaseq_proteome_matchedCohort_featureCounts_supTempGyrus.rds')
```

``` r
library(edgeR)
y <- DGEList(counts = mycounts$counts, genes = mycounts$annotation)
keep <- filterByExpr(y)
table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ##  8403 17299

``` r
y <- y[keep,,keep.lib.sizes = FALSE]
y <- normLibSizes(y)
rna_rpkm <- rpkm(y) |> as.data.frame() |> 
    rownames_to_column() |> 
    mutate(rowname = as.numeric(rowname)) |> 
    left_join(annotables::grch37, by = c('rowname' = 'entrez'))
```

``` r
library(ggpubr)
library(ggpmisc)
rnaoi <- rna_rpkm %>% 
    dplyr::filter(symbol %in% c('MAPT', 'SNCA')) %>%
    dplyr::select(symbol, ends_with('.bam')) %>% 
    unique() %>% 
    pivot_longer(-symbol) %>% 
    pivot_wider(names_from = 'symbol') %>% 
    left_join(metad_rna, by = 'name') %>% 
    left_join(metad, by = c('individualID' = 'individualIdentifier')) %>% 
    unique()
rnaoi %>% 
    # dplyr::filter(MAPT > 20) %>%
    ggplot(aes(x = MAPT, y = SNCA)) +
    geom_point(aes(color = as.factor(bbscore))) +
    geom_smooth(method = 'lm') +
    stat_cor(label.y.npc = 'bottom') +
    scale_color_brewer(guide = 'none', palette = 'YlOrRd') +
    labs(title = 'MSBB RNA-seq', y = 'SNCA [rpkm]', x = 'MAPT [rpkm]')
```

![](gh_04_MSBB_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Get genotypes from WGS:

``` r
library(VariantAnnotation)
metad_wgs <- read_tsv('/earth/public_data/datasets/MSBB/Metadata/AMP-AD_MSBB_WGS__sample_barcode_brainBankID..tsv') %>% left_join(metad, by = 'individualIdentifier') %>% 
    mutate(sampleIdentifier = as.character(sampleIdentifier))
rs356182 <- as_tibble(geno(readVcf('/earth/public_data/datasets/MSBB/genomicVariants/MSBB_rs356182_grch37.vcf'))$GT) %>% 
    pivot_longer(everything()) %>% 
    mutate(allele = case_when(value %in% c('1|0', '0|1') ~ 'HET',
                              value == '0|0' ~ 'RISK',
                              value == '1|1' ~ 'PROT')) %>% 
    left_join(metad_wgs,  by = c('name' = 'sampleIdentifier')) %>% 
    left_join(proteome_cor, by = c('name.y' = 'name')) %>% 
    dplyr::filter(!is.na(SNCA))
rs356168 <- as_tibble(geno(readVcf('/earth/public_data/datasets/MSBB/genomicVariants/MSBB_rs356168_grch37.vcf'))$GT) %>% 
    pivot_longer(everything()) %>% 
    mutate(allele = case_when(value %in% c('1|0', '0|1') ~ 'HET',
                              value == '0|0' ~ 'RISK',
                              value == '1|1' ~ 'PROT')) %>% 
    left_join(metad_wgs,  by = c('name' = 'sampleIdentifier')) %>% 
    left_join(proteome_cor, by = c('name.y' = 'name')) %>% 
    dplyr::filter(!is.na(SNCA))
rowRanges(readVcf('/earth/public_data/datasets/MSBB/genomicVariants/MSBB_rs356182_grch37.vcf')) # REF is G, ALT is A
```

    ## GRanges object with 1 range and 5 metadata columns:
    ##                  seqnames    ranges strand | paramRangeID            REF
    ##                     <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet>
    ##   4:90626111_G/A        4  90626111      * |           NA              G
    ##                                 ALT      QUAL      FILTER
    ##                  <DNAStringSetList> <numeric> <character>
    ##   4:90626111_G/A                  A        NA           .
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
rowRanges(readVcf('/earth/public_data/datasets/MSBB/genomicVariants/MSBB_rs356168_grch37.vcf')) # REF is G, ALT is A
```

    ## GRanges object with 1 range and 5 metadata columns:
    ##                  seqnames    ranges strand | paramRangeID            REF
    ##                     <Rle> <IRanges>  <Rle> |     <factor> <DNAStringSet>
    ##   4:90674431_G/A        4  90674431      * |           NA              G
    ##                                 ALT      QUAL      FILTER
    ##                  <DNAStringSetList> <numeric> <character>
    ##   4:90674431_G/A                  A        NA           .
    ##   -------
    ##   seqinfo: 1 sequence from an unspecified genome; no seqlengths

``` r
# A is protective for both
```

``` r
rnaoi %>% left_join(rs356168, by = c('individualID' =  'individualIdentifier.x') ) %>% 
    dplyr::select(allele, MAPT.x, SNCA.x) %>% 
    dplyr::filter(!is.na(allele)) %>% 
    mutate(allele = factor(allele, levels = c('PROT', 'HET', 'RISK'))) %>% 
    ggplot(aes(x = allele, y = SNCA.x, fill = allele)) +
    geom_violin() +
    geom_jitter(width = 0.01, alpha = 0.5) +
    stat_compare_means(method = 'wilcox', method.args = list(alternative = 'greater'), label = 'p.signif', size = 6, comparisons = list(c('HET', 'PROT'), c('HET', 'RISK'))) +
    # stat_compare_means(method = 't.test', method.args = list(alternative = 'two.sided'), label = 'p.label', size = 6, comparisons = list(c('PROT_HET', 'RISK'))) +
    scale_fill_discrete(guide = FALSE) +
    ylim(5, 85) +
    labs(title = "Allele-specific SNCA expression", x = '', y = 'SNCA [rpkm]', fill = '')
```

![](gh_04_MSBB_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## 0.3 Session info

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Arch Linux
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/liblapack.so.3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] VariantAnnotation_1.52.0    Rsamtools_2.22.0           
    ##  [3] Biostrings_2.74.1           XVector_0.46.0             
    ##  [5] SummarizedExperiment_1.36.0 Biobase_2.66.0             
    ##  [7] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
    ##  [9] IRanges_2.40.1              S4Vectors_0.44.0           
    ## [11] MatrixGenerics_1.18.1       matrixStats_1.5.0          
    ## [13] BiocGenerics_0.52.0         ggpmisc_0.6.1              
    ## [15] ggpp_0.5.8-1                edgeR_4.4.2                
    ## [17] limma_3.62.2                Rsubread_2.20.0            
    ## [19] ggpubr_0.6.0                httr_1.4.7                 
    ## [21] viridis_0.6.5               viridisLite_0.4.2          
    ## [23] RColorBrewer_1.1-3          patchwork_1.3.0            
    ## [25] ggraph_2.2.1                lubridate_1.9.4            
    ## [27] forcats_1.0.0               stringr_1.5.1              
    ## [29] dplyr_1.1.4                 purrr_1.0.4                
    ## [31] readr_2.1.5                 tidyr_1.3.1                
    ## [33] tibble_3.2.1                ggplot2_3.5.1              
    ## [35] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.17.1        jsonlite_1.8.9           magrittr_2.0.3          
    ##  [4] GenomicFeatures_1.58.0   farver_2.1.2             rmarkdown_2.29          
    ##  [7] BiocIO_1.16.0            zlibbioc_1.52.0          vctrs_0.6.5             
    ## [10] memoise_2.0.1            RCurl_1.98-1.16          rstatix_0.7.2           
    ## [13] htmltools_0.5.8.1        S4Arrays_1.6.0           polynom_1.4-1           
    ## [16] curl_6.2.0               broom_1.0.7              SparseArray_1.6.1       
    ## [19] Formula_1.2-5            cachem_1.1.0             GenomicAlignments_1.42.0
    ## [22] igraph_2.1.4             lifecycle_1.0.4          pkgconfig_2.0.3         
    ## [25] Matrix_1.7-2             R6_2.6.0                 fastmap_1.2.0           
    ## [28] GenomeInfoDbData_1.2.13  digest_0.6.37            colorspace_2.1-1        
    ## [31] AnnotationDbi_1.68.0     RSQLite_2.3.9            labeling_0.4.3          
    ## [34] timechange_0.3.0         polyclip_1.10-7          abind_1.4-8             
    ## [37] mgcv_1.9-1               compiler_4.4.2           bit64_4.6.0-1           
    ## [40] withr_3.0.2              backports_1.5.0          BiocParallel_1.38.0     
    ## [43] carData_3.0-5            DBI_1.2.3                ggforce_0.4.2           
    ## [46] ggsignif_0.6.4           MASS_7.3-64              quantreg_6.00           
    ## [49] DelayedArray_0.32.0      rjson_0.2.23             tools_4.4.2             
    ## [52] glue_1.8.0               restfulr_0.0.15          nlme_3.1-167            
    ## [55] grid_4.4.2               generics_0.1.3           BSgenome_1.74.0         
    ## [58] gtable_0.3.6             tzdb_0.4.0               hms_1.1.3               
    ## [61] tidygraph_1.3.1          car_3.1-3                utf8_1.2.4              
    ## [64] ggrepel_0.9.6            pillar_1.10.1            vroom_1.6.5             
    ## [67] splines_4.4.2            tweenr_2.0.3             lattice_0.22-6          
    ## [70] rtracklayer_1.66.0       survival_3.8-3           bit_4.5.0.1             
    ## [73] SparseM_1.84-2           tidyselect_1.2.1         locfit_1.5-9.11         
    ## [76] knitr_1.49               annotables_0.2.0         gridExtra_2.3           
    ## [79] xfun_0.50                graphlayouts_1.2.2       statmod_1.5.0           
    ## [82] stringi_1.8.4            UCSC.utils_1.2.0         yaml_2.3.10             
    ## [85] evaluate_1.0.3           codetools_0.2-20         cli_3.6.4               
    ## [88] munsell_0.5.1            Rcpp_1.0.14              png_0.1-8               
    ## [91] XML_3.99-0.18            parallel_4.4.2           MatrixModels_0.5-3      
    ## [94] blob_1.2.4               bitops_1.0-9             scales_1.3.0            
    ## [97] crayon_1.5.3             rlang_1.1.5              KEGGREST_1.46.0
