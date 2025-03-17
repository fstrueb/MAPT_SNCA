WGS analyses pertaining to Figure 1
================
<felix.struebing@med.uni-muenchen.de>
2025-03-17

- [1 Prerequisites](#1-prerequisites)
  - [1.1 Get patient data](#11-get-patient-data)
  - [1.2 Filtering](#12-filtering)
  - [1.3 N of samples](#13-n-of-samples)
- [2 Analysis](#2-analysis)
  - [2.1 Fig. 1b](#21-fig-1b)
    - [2.1.1 Assembly](#211-assembly)
  - [2.2 PRSKB](#22-prskb)
    - [2.2.1 Fig. 1c](#221-fig-1c)
    - [2.2.2 Fig. 1d](#222-fig-1d)
    - [2.2.3 Supp. Fig. 1](#223-supp-fig-1)
    - [2.2.4 Fig. 1e](#224-fig-1e)
    - [2.2.5 Fig. 1f](#225-fig-1f)
  - [2.3 GBA mutation testing](#23-gba-mutation-testing)
  - [2.4 Familial AD mutations](#24-familial-ad-mutations)
    - [2.4.1 Fig. 1g](#241-fig-1g)
  - [2.5 Assembly](#25-assembly)
- [3 Supp. Table 1](#3-supp-table-1)
  - [3.1 Gather PRS](#31-gather-prs)
    - [3.1.1 AD](#311-ad)
    - [3.1.2 DLB](#312-dlb)
    - [3.1.3 PD](#313-pd)
  - [3.2 Gather MAPT haplotype](#32-gather-mapt-haplotype)
  - [3.3 Gather GBA mutations](#33-gather-gba-mutations)
  - [3.4 Gather SNCA alleles](#34-gather-snca-alleles)
- [4 Session Info](#4-session-info)

``` r
library(tidyverse)
library(broom)
library(ggpubr)
library(ggsci)
library(scales)
library(patchwork)
library(VariantAnnotation)
AD_color <- pal_igv('alternating')(2)[1]
ASYN_color <- pal_igv('alternating')(2)[2]
WLDY_color <- '#66cc00'
source('ggplot_theme_patchwork.R')
source('scripts/wgs_functions.R')
```

# 1 Prerequisites

## 1.1 Get patient data

``` r
rzids  <- read_tsv('resources/synergy_cohort_metad.tsv') %>% 
    mutate(asyn_bin = factor(ifelse(asyn == 'yes', 1, 0), levels = c(0, 1), labels = c('no', 'yes')))
```

## 1.2 Filtering

Remove familial AD mutations, brains with BB \< 4, and one sample with
bad quality metrics:

``` r
excludesamples  <- rzids  %>% 
    filter(!is.na(variant) | braakbraak <= 3 | sample == 'Sample_57')  %>% 
    dplyr::pull(sample)
```

## 1.3 N of samples

The superior frontal cortex was stained with antibodies against asyn
(clone 42). Positivity, equivalent to Braak LBP stage 6 (McKeith
neocortical), is coded as “yes” or “no” in the “asyn” variable.

``` r
# number of samples
rzids %>% filter(!sample %in% excludesamples)  %>% nrow(.)
```

    ## [1] 135

``` r
# distribution of samples
rzids %>% filter(!sample %in% excludesamples) %>% 
    dplyr::count(asyn)
```

    ## # A tibble: 2 × 2
    ##   asyn      n
    ##   <chr> <int>
    ## 1 no      100
    ## 2 yes      35

# 2 Analysis

## 2.1 Fig. 1b

``` r
p1b1  <- rzids  %>% 
    filter(!sample %in% excludesamples) %>% 
    ggplot(aes(x = braakbraak, fill = asyn)) +
    geom_bar(position = 'dodge', color = 'black') +
    coord_flip() +
    scale_fill_manual(values = c(AD_color, ASYN_color)) +
    labs(title = "Braak tau stages", x = 'Braak tau stage', y = 'n cases', fill = 'aSyn') +
    ylim(0, 55) +
    theme_patch() +
    theme(legend.position = 'top')
```

``` r
p1b2 <- rzids  %>% 
    filter(!sample %in% excludesamples) %>% 
    group_by(asyn, APOE_WGS) %>% 
    tally() %>% 
    ggplot(aes(x = asyn, y = n, fill = APOE_WGS)) +
    geom_col(position = 'dodge', color = 'black') +
    stat_compare_means(comparisons = list(c('no', 'yes')), method = 'wilcox.test', label = 'p.signif') +
    labs(title = "WGS cohort: APOE genotypes", x = '', y = 'n cases', fill = 'APOE') +
    ylim(0, 42) +
    theme_patch()
```

### 2.1.1 Assembly

``` r
p1b <- p1b1 | p1b2
# pdf('Fig_1b.pdf', height = 4, width = 6)
p1b
```

![](v2_Figure_1_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# dev.off()
```

## 2.2 PRSKB

The bash command used was
`sh /opt/PrskbCLITool/runPrsCLI.sh -f wellderly200_synergy_wgs_merged_hg38_uniquevars_rsidanno.vcf -o prskb_wellderly_synergy_output.tsv -r hg38 -c 0.05 -p EUR`
after annotating the joint (`bcftools merge`) vcf with
`bcftools annotate --threads 24 -a /earth/public_data/vcf/dbSNP_151_GRCh38_GATK.vcf.gz -c ID -Oz -o wellderly200_synergy_wgs_merged_hg38_uniquevars_rsidanno.vcf wellderly200_synergy_wgs_merged_hg38_uniquevars_noalt.vcf.gz`.

``` r
prskb <- read_tsv('/neptune/wellderly/vcf/prskb_wellderly_synergy_output.tsv')
```

Data wrangling:

``` r
prskb_long <- prskb %>% 
    pivot_longer(cols = c(starts_with('Sample'), starts_with('23'))) %>% 
    left_join(rzids, by = c('name' = 'sample')) %>% 
    mutate(group = case_when(
        asyn == 'no' ~ 'AD',
        asyn == 'yes' ~ 'AD+ASYN',
        is.na(asyn) ~ 'WLDY'
    )) %>% 
    dplyr::filter(!name %in% excludesamples) %>% 
    mutate(asyn_bin = factor(ifelse(asyn == 'yes', 1, 0), levels = c(0, 1), labels = c('no', 'yes'))) |> 
    filter(value != 'NF') |> 
    mutate(value_num = as.numeric(value)) |> 
    mutate(study_trait = paste(`Study ID`, `Reported Trait`, Trait, `P-Value Annotation`, sep = '_')) |> 
    dplyr::select(study_trait, name, value_num, APOE_WGS, sex, group, asyn_bin) |> 
    filter(is.na(APOE_WGS) | APOE_WGS != 'E2/E4') |> #otherwise unequal factor levels
    mutate(APOE_WGS = droplevels(as.factor(APOE_WGS))) 
# missing: 02-21-25 necessary for sex correction
```

### 2.2.1 Fig. 1c

``` r
prs_ad <- prskb_long |> 
    # filter(study_trait %in% prskb_sig$study_trait) |> 
    # dplyr::filter(grepl('Park', study_trait)) |> 
    dplyr::filter(study_trait == "GCST007320_Alzheimer's Disease Or Family History Of Alzheimer's Disease_Alzheimer Disease_NA") %>% 
    dplyr::mutate(asyn_bin = ifelse(asyn_bin == 'yes', 'AD+ASYN', 'AD')) |> 
    mutate(beta_zscore = scale(value_num)) |>
    mutate(group = factor(group, levels = c('WLDY', 'AD', 'AD+ASYN'))) %>% 
    ggplot(aes(x = group, y = beta_zscore, fill = group)) +
    geom_violin() +
    geom_jitter(width = 0.01, alpha= 0.5) +
    stat_compare_means(label = 'p.signif', comparisons = list(
        c('WLDY', 'AD'), c('WLDY', 'AD+ASYN'), c('AD', 'AD+ASYN')), 
        size = 6) +
    scale_fill_manual(values = c(WLDY_color, AD_color, ASYN_color), guide = FALSE) +
    ylim(-3, 5.5) +
    theme_patch() +
    labs(title = "PRS: Late-onset AD", x = '', y = 'beta [z-scored]', subtitle = 'Jansen et al. 2019', fill = '')
prs_ad
```

![](v2_Figure_1_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### 2.2.2 Fig. 1d

``` r
prs_pd  <- prskb_long |> 
    dplyr::filter(!name %in% excludesamples) |> 
    dplyr::filter(study_trait == "GCST009325_Parkinson's Disease Or First Degree Relation To Individual With Parkinson's Disease_Parkinson Disease_NA") |> 
    # dplyr::mutate(asyn_bin = ifelse(asyn_bin == 'yes', 'AD+ASYN', 'AD')) |> 
    mutate(beta_zscore = scale(value_num)) |> 
    mutate(group = factor(group, levels = c('WLDY', 'AD', 'AD+ASYN'))) %>% 
    ggplot(aes(x = group, y = beta_zscore, fill = group)) +
    geom_violin() +
    geom_jitter(width = 0.01, alpha= 0.5) +
    stat_compare_means(label = 'p.signif', comparisons = list(
        c('WLDY', 'AD'), c('WLDY', 'AD+ASYN'), c('AD', 'AD+ASYN')), 
        size = 6) +
    scale_fill_manual(values = c(WLDY_color, AD_color, ASYN_color), guide = FALSE) +
    ylim(-3, 4.5) +
    theme_patch() +
    labs(title = "PRS: Parkinson's disease", x = '', y = 'beta [z-scored]', subtitle = 'Nalls et al. 2019', fill = '')
prs_pd
```

![](v2_Figure_1_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### 2.2.3 Supp. Fig. 1

PRS for DLB:

``` r
prs_dlb <- prskb_long |> 
    dplyr::filter(!name %in% excludesamples) |> 
    dplyr::filter(study_trait == "GCST90001390_Dementia With Lewy Bodies_Lewy Body Dementia_NA") |> 
    # dplyr::mutate(asyn_bin = ifelse(asyn_bin == 'yes', 'AD+ASYN', 'AD')) |> 
    mutate(beta_zscore = scale(value_num)) |> 
    mutate(group = factor(group, levels = c('WLDY', 'AD', 'AD+ASYN'))) %>% 
    ggplot(aes(x = group, y = beta_zscore, fill = group)) +
    geom_violin() +
    geom_jitter(width = 0.01, alpha= 0.5) +
    stat_compare_means(label = 'p.signif', comparisons = list(
        c('WLDY', 'AD'), c('WLDY', 'AD+ASYN'), c('AD', 'AD+ASYN')), 
        size = 6) +
    scale_fill_manual(values = c(WLDY_color, AD_color, ASYN_color), guide = FALSE) +
    ylim(-2, 7) +
    theme_patch() +
    labs(title = "PRS: Dementia with Lewy Bodies (DLB)", x = '', y = 'beta [z-scored]', subtitle = 'Chia et al. 2021', fill = '')
# prs_dlb
```

``` r
# pdf('SFig1_DLB_PRS.pdf', width = 6, height = 4)
prs_dlb
```

![](v2_Figure_1_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# dev.off()
```

### 2.2.4 Fig. 1e

``` r
#### we used GCST009325_Parkinson's Disease Or First Degree Relation To Individual With Parkinson's Disease_Parkinson Disease_NA
prskb_full <- read_tsv('/earth/AD_WGS_backup/complete_genotypes_VQSR/prskb_output_full.tsv') 
prskb_full_test <- prskb_full %>% filter(`Study ID` == 'GCST009325') %>% 
    mutate(prs = as.numeric(`Polygenic Risk Score`)) %>% 
    arrange(-prs) %>% 
    left_join(rzids, by = c('Sample' = 'sample')) %>% 
    dplyr::select(Sample, prs, asyn, Percentile, protvars = `Protective Variants`, riskvars = `Risk Variants`, rzid, age, APOE_WGS, braakbraak, sex) %>% 
    filter(!Sample %in% excludesamples) %>% 
    left_join(., read_delim('../AD_WGS/MAPT_haplotypes.tsv', delim = '\t') %>% dplyr::select(name, MAPT_haplotype = `case_when(...)`), by = c('Sample' = 'name'))
```

``` r
regr_on_vars <- prskb_full_test %>% 
    mutate(rv = strsplit(riskvars, "\\|")) %>%
    unnest(rv) %>%
    group_by(Sample) %>%
    mutate(row = paste0('risk_', row_number())) %>%
    spread(row, rv) %>% 
    mutate(pv = strsplit(protvars, "\\|")) %>%
    unnest(pv) %>%
    group_by(Sample) %>%
    mutate(row = paste0('prot_', row_number())) %>%
    spread(row, pv) %>% 
    dplyr::select(Sample, prs, asyn, Percentile, starts_with('risk_'), starts_with('prot_')) %>% 
    pivot_longer(cols = starts_with('risk_'), names_to = 'riskvars', values_to = 'risk_rsid') %>% 
    pivot_longer(cols = starts_with('prot_'), names_to = 'protvars', values_to = 'prot_rsid') %>% 
    dplyr::select(-riskvars, -protvars) %>% 
    pivot_longer(cols = c(risk_rsid, prot_rsid)) %>% 
    unique() %>% 
    na.omit() %>% 
    mutate(Percentile = ifelse(Percentile == '74-75', '75', Percentile)) %>% 
    mutate(Percentile = as.numeric(Percentile))
```

Make sure there are no intersections between risk and protective
variants:

``` r
temp_prot_vars <- regr_on_vars %>% filter(name == 'prot_rsid') %>% dplyr::pull(value) %>% unique()
temp_risk_vars <- regr_on_vars %>% filter(name == 'risk_rsid') %>% dplyr::pull(value) %>% unique()
intersect(temp_prot_vars, temp_risk_vars)
```

    ## character(0)

Fit a linear model: PRS percentile on number of risk/protective
variants:

``` r
regr_on_vars_res <- tidy(summary(lm(Percentile ~ value, data = regr_on_vars))) %>% 
    filter(term != '(Intercept)') %>% 
    arrange(p.value) %>% 
    mutate(rsid = str_sub(term, start = 6L)) %>% 
    mutate(var_mod = ifelse(rsid %in% temp_prot_vars, 'protective', 'risk'))
# regr_on_vars_res %>% filter(rsid == 'rs356182')
```

``` r
p1e <- regr_on_vars_res %>% 
    filter(p.value < 0.05) %>% 
    mutate(locus = ifelse(var_mod == 'protective', '17q21.31', '1q22')) %>% 
    dplyr::select(rsid, estimate, p.value, var_mod, locus, std.error) %>% 
    ggplot(aes(x = reorder(rsid, -(estimate)), y = estimate/100, fill = locus, ymin = estimate/100 - std.error/100, ymax = estimate/100 + std.error/100)) +
    geom_col(color = 'black') +
    geom_errorbar(width = 0.2) +
    labs(x = 'Variant', y = 'Estimate (PRS percentile ~ variants)', title = 'Overrepresented variants in AD+ASYN') +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() +
    theme_patch() +
    theme(legend.position = 'right')
```

### 2.2.5 Fig. 1f

``` r
temp_haplo <- prskb_full_test %>% 
    dplyr::select(asyn, prs, MAPT_haplotype) %>% 
    mutate(asyn_bin = ifelse(asyn == 'yes', 1, 0))
tidy(chisq.test(x = temp_haplo$asyn_bin, y = temp_haplo$MAPT_haplotype))
```

    ## # A tibble: 1 × 4
    ##   statistic p.value parameter method                    
    ##       <dbl>   <dbl>     <int> <chr>                     
    ## 1      8.42  0.0149         2 Pearson's Chi-squared test

``` r
# tidy(kruskal.test(asyn_bin ~ MAPT_haplotype, data = temp_haplo))
# tidy(glm(asyn_bin ~ 0 + MAPT_haplotype, data = temp_haplo, family = binomial('logit')))
```

``` r
p1f <- temp_haplo %>%
    mutate(myhaplo = ifelse(MAPT_haplotype == 'H1', 'H1', 'HET/H2')) %>%
    # filter(myhaplo == 'H1') %>% 
    mutate(beta_zscore = scale(prs)) |>
    mutate(group = ifelse(asyn == 'yes', 'AD+ASYN', 'AD')) %>%
    mutate(group = factor(group, levels = c('AD', 'AD+ASYN'))) %>%
    group_by(group, myhaplo) %>% 
    tally() %>% 
    ggplot(aes(x = group, y = n, fill = myhaplo)) +
    geom_col(position = 'fill', color = 'black') +
    scale_fill_manual(values = c('grey50', 'grey90')) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = '', y = '', title = 'MAPT haplotypes per group') +
    theme_patch()
```

<!-- ```{r} -->
<!-- prskb_full_test %>%  -->
<!--   mutate(rv = strsplit(riskvars, "\\|")) %>% -->
<!--         unnest(rv) %>% -->
<!--         group_by(Sample) %>% -->
<!--         mutate(row = paste0('risk_', row_number())) %>% -->
<!--         spread(row, rv) %>%  -->
<!--   mutate(pv = strsplit(protvars, "\\|")) %>% -->
<!--         unnest(pv) %>% -->
<!--         group_by(Sample) %>% -->
<!--         mutate(row = paste0('prot_', row_number())) %>% -->
<!--         spread(row, pv) %>%  -->
<!--   dplyr::select(Sample, prs, asyn, Percentile, starts_with('risk_'), starts_with('prot_')) %>%  -->
<!--   pivot_longer(cols = starts_with('risk_'), names_to = 'riskvars', values_to = 'risk_rsid') %>%  -->
<!--   pivot_longer(cols = starts_with('prot_'), names_to = 'protvars', values_to = 'prot_rsid') %>%  -->
<!--   dplyr::select(-riskvars, -protvars) %>%  -->
<!--   pivot_longer(cols = c(risk_rsid, prot_rsid)) %>%  -->
<!--   unique() %>%  -->
<!--   na.omit() %>%  -->
<!--   group_by(asyn, name) %>%  -->
<!--   tally() %>%  -->
<!--   pivot_wider(values_from = n) %>%  -->
<!--   left_join(rzids, by = 'asyn') %>%  -->
<!--   filter(!sample %in% excludesamples) %>%  -->
<!--   group_by(asyn) %>%  -->
<!--   add_tally(name = 'n_in_group') %>%  -->
<!--   dplyr::select(asyn, n_in_group, prot_rsid, risk_rsid) %>%  -->
<!--   unique() %>%  -->
<!--   mutate(norm_prot = prot_rsid / n_in_group, -->
<!--          norm_risk = risk_rsid / n_in_group) -->
<!-- ``` -->

## 2.3 GBA mutation testing

``` r
# metad <- readxl::read_xlsx('resources/AD_Fälle_NBM_FLS.xlsx', sheet = 1)
# alzmut <- readRDS('resources/alzforum_mutations_scraped_220921.rds') %>% 
#     mutate(mutation_alias = str_extract(Mutation, '\\(+.*')) %>% 
#     mutate(Mutation = str_remove(Mutation, '\\(+.*'))
myvcf <- '/earth/AD_WGS_backup/complete_genotypes_VQSR/ADWGS_complete_recal_indels_snps.vcf.gz' 
```

## 2.4 Familial AD mutations

``` r
mut_list <- invisible(lapply(c('GBA'), function(x) {extract_genotypes(x, myvcf, 'chm13') %>% 
        dplyr::mutate(query_gene = x)}))
```

    ## [1] "GBA"

``` r
vcf_mutations <- rbindlist(mut_list) %>% as_tibble()
```

Annotate alleles with genotyping quality for hard filtering (can skip
after VQSR) and select MANE isoform:

``` r
temp_myres <- vcf_mutations %>% 
    dplyr::select(rowname) %>% 
    mutate(pos = str_remove(rowname, pattern = '[\\_]+.*')) %>% 
    dplyr::pull(pos) %>% unique()
quals <- lapply(temp_myres, function(x) {get_vcfpos_quals(myvcf, x)}) %>% rbindlist(.)
myres <- vcf_mutations %>% 
    # dplyr::filter(consequence != 'synonymous') %>% 
    mutate(pos = str_remove(rowname, pattern = '[\\_]+.*')) %>% 
    left_join(., quals, by = c('pos' = 'position')) %>% 
    # dplyr::filter(decision == 'OK') %>% 
    dplyr::filter(grepl('MANE', gencode_tag))
```

Discard synonymous mutations:

``` r
myres_clean <- myres %>% 
    # dplyr::filter(decision == 'OK') %>% 
    dplyr::filter(consequence != 'synonymous') %>% 
    mutate(mutation_aa = paste0(aaref, aapos, aaalt)) %>% 
    mutate(mut_alias = paste0('(', str_replace(str_sub(str_extract(rowname, pattern = '[\\_]+.*'), start = 2L), pattern = '\\/', replacement = '>'), ')')) %>% 
    dplyr::select(sample = name, hs1_position_allele = rowname, allele = value, consequence, gene_symbol = query_gene, mutation_protein = mutation_aa, mutation_gene = mut_alias) %>% 
    left_join(., rzids %>% dplyr::select(sample, asyn) , by = 'sample')
knitr::kable(myres_clean)
```

| sample     | hs1_position_allele | allele | consequence   | gene_symbol | mutation_protein | mutation_gene | asyn |
|:-----------|:--------------------|:-------|:--------------|:------------|:-----------------|:--------------|:-----|
| Sample_24  | chr1:154373809_A/G  | 0\|1   | nonsynonymous | GBA         | L483P            | (A\>G)        | yes  |
| Sample_20  | chr1:154374325_C/T  | 0\|1   | nonsynonymous | GBA         | R434H            | (C\>T)        | no   |
| Sample_155 | chr1:154374400_T/C  | 1\|0   | nonsynonymous | GBA         | N409S            | (T\>C)        | yes  |
| Sample_34  | chr1:154374803_G/A  | 0\|1   | nonsynonymous | GBA         | T408M            | (G\>A)        | no   |
| Sample_50  | chr1:154374803_G/A  | 0\|1   | nonsynonymous | GBA         | T408M            | (G\>A)        | no   |
| Sample_82  | chr1:154374803_G/A  | 0\|1   | nonsynonymous | GBA         | T408M            | (G\>A)        | no   |
| Sample_90  | chr1:154374803_G/A  | 0\|1   | nonsynonymous | GBA         | T408M            | (G\>A)        | yes  |
| Sample_159 | chr1:154374933_C/T  | 1\|0   | nonsynonymous | GBA         | E365K            | (C\>T)        | yes  |
| Sample_29  | chr1:154374933_C/T  | 1\|0   | nonsynonymous | GBA         | E365K            | (C\>T)        | no   |
| Sample_31  | chr1:154374933_C/T  | 0\|1   | nonsynonymous | GBA         | E365K            | (C\>T)        | no   |
| Sample_39  | chr1:154374933_C/T  | 1\|0   | nonsynonymous | GBA         | E365K            | (C\>T)        | yes  |
| Sample_44  | chr1:154374933_C/T  | 0\|1   | nonsynonymous | GBA         | E365K            | (C\>T)        | no   |
| Sample_65  | chr1:154374933_C/T  | 1\|0   | nonsynonymous | GBA         | E365K            | (C\>T)        | no   |
| Sample_67  | chr1:154374933_C/T  | 0\|1   | nonsynonymous | GBA         | E365K            | (C\>T)        | no   |
| Sample_162 | chr1:154376726_G/T  | 0\|1   | nonsynonymous | GBA         | D242E            | (G\>T)        | yes  |

### 2.4.1 Fig. 1g

rs356182: REF (G) is risk. rs356168: ALT (A) is protective.

``` r
myvcf <- readVcf('/earth/AD_WGS_backup/complete_genotypes_VQSR/ADWGS_rsid_hg38_snca_elements_eqtls_v2.vcf')
mygenos_fil <- geno(myvcf)$GT %>% 
    as.data.frame(.) %>% 
    rownames_to_column(., 'rsid') %>% 
    mutate(alt_allele = as.character(unlist(alt(myvcf)))) %>% 
    mutate(ref_allele = as.character(ref(myvcf))) %>% 
    pivot_longer(c(-rsid, -alt_allele, -ref_allele)) %>% 
    right_join(., rzids, by = c('name' = 'sample')) %>% 
    dplyr::filter(rsid %in% c('rs356182', 'rs356168')) %>% 
    mutate(risk_allele = case_when(
        rsid == 'rs356182' & value == '1|1' ~ 'protective',
        rsid == 'rs356182' & value == '0|0' ~ 'risk',
        rsid == 'rs356168' & value == '1|1' ~ 'protective',
        rsid == 'rs356168' & value == '0|0' ~ 'risk',
        value == '1|0' | value == '0|1' ~ 'HET'
    )) %>% 
    dplyr::filter(!name %in% excludesamples) %>% 
    mutate(asyn_bin = as.factor(asyn)) %>% 
    mutate(risk_allele = as.factor(risk_allele))
mygenos_fil$risk_allele <- factor(mygenos_fil$risk_allele, levels = c('protective', 'HET', 'risk'))
```

Fisher’s test after concatenating HET and protective:

``` r
mygenos_fish <- mygenos_fil %>% 
    mutate(asyn = ifelse(asyn == 'yes', 'asyn', 'no')) %>% 
    mutate(risk_group = ifelse(risk_allele == 'risk', 'high_risk', 'low_risk')) %>% 
    group_by(rsid, asyn, risk_group) %>% 
    tally() %>% 
    split(., f = .$rsid) %>% 
    lapply(., function(x) {
        x %>% ungroup() %>% 
            dplyr::select(asyn, risk_group, n) %>% 
            pivot_wider(id_cols = asyn, names_from = risk_group, values_from = n) %>% 
            column_to_rownames('asyn') %>% 
            as.matrix()
    })
mygenos_fish
```

    ## $rs356168
    ##      high_risk low_risk
    ## asyn         9       26
    ## no          10       90
    ## 
    ## $rs356182
    ##      high_risk low_risk
    ## asyn         6       29
    ## no           6       94

``` r
mygenos_tested <- lapply(mygenos_fish, fisher.test, alternative = 'two.sided')
```

``` r
mygenos_tested_plot <- lapply(mygenos_tested, function(x) {
    ci <- x$conf.int
    or <- x$estimate
    df <- data.frame(ci_min = ci[1], ci_max = ci[2], or = or)
}) %>% data.table::rbindlist(idcol = 'SNP') 
```

``` r
p1g <- mygenos_tested_plot %>% 
    ggplot(aes(x = or, xmax = ci_max, xmin = ci_min, y = SNP)) +
    geom_errorbarh(size = .5, height = .2, color = 'gray50') +
    geom_point(size = 3.5, color = 'gray20') +
    geom_vline(xintercept = 1, size = .25, linetype = 'dashed') +
    labs(x = 'OR [log10]', y = 'Variant', title = 'Odds Ratios: SNCA variants') +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
    xlim(c(0.1, 20)) +
    # xlim(-20, 20) +
    theme_patch()
```

## 2.5 Assembly

``` r
p1cde_layout <- "
##AABB
CCDDEE
CCDDFF
CCDDGG
"
p1cde <- wrap_plots(p1b1, p1b2, prs_ad, prs_pd, p1e, p1f, p1g, design = p1cde_layout)
# p1cde
```

``` r
# pdf('Fig1.pdf', width = 14, height = 12)
# p1cde
# dev.off()
```

# 3 Supp. Table 1

## 3.1 Gather PRS

### 3.1.1 AD

``` r
stab1_ad <- prskb_long %>% 
    dplyr::filter(study_trait == "GCST007320_Alzheimer's Disease Or Family History Of Alzheimer's Disease_Alzheimer Disease_NA") %>% 
    mutate(beta_zscore = scale(value_num)) %>% 
    filter(group %in% c('AD', 'AD+ASYN')) %>% 
    dplyr::select(sample = name, AD_PRS = beta_zscore) %>% 
    mutate(AD_PRS = as.numeric(AD_PRS))
```

### 3.1.2 DLB

``` r
stab1_dlb <- prskb_long %>% 
    dplyr::filter(study_trait == "GCST90001390_Dementia With Lewy Bodies_Lewy Body Dementia_NA") %>% 
    mutate(beta_zscore = scale(value_num))  %>% 
    filter(group %in% c('AD', 'AD+ASYN')) %>% 
    dplyr::select(sample = name, DLB_PRS = beta_zscore) %>% 
    mutate(DLB_PRS = as.numeric(DLB_PRS))
```

### 3.1.3 PD

``` r
stab1_pd <- prskb_long %>%  
    dplyr::filter(study_trait == "GCST009325_Parkinson's Disease Or First Degree Relation To Individual With Parkinson's Disease_Parkinson Disease_NA") %>% 
    mutate(beta_zscore = scale(value_num)) %>% 
    filter(group %in% c('AD', 'AD+ASYN')) %>% 
    dplyr::select(sample = name, PD_PRS = beta_zscore) %>% 
    mutate(PD_PRS = as.numeric(PD_PRS))
```

## 3.2 Gather MAPT haplotype

``` r
stab1_mapt <- prskb_full_test  %>% 
    dplyr::select(sample = Sample, MAPT_haplotype) %>% 
    mutate(MAPT_haplotype = ifelse(MAPT_haplotype == 'H2_het', 'HET', MAPT_haplotype)) %>% 
    mutate(MAPT_haplotype = ifelse(MAPT_haplotype == 'H2_hom', 'H2', MAPT_haplotype)) 
```

## 3.3 Gather GBA mutations

``` r
stab1_gba <- myres_clean %>% 
    dplyr::select(sample, GBA1_mutation = mutation_protein)
```

## 3.4 Gather SNCA alleles

``` r
stab1_snca <- mygenos_fil %>% 
    dplyr::select(sample = name, rsid, risk_allele) %>% 
    pivot_wider(names_from = rsid, values_from = risk_allele)
```

``` r
rzids_export <- rzids %>% 
    filter(!sample %in% excludesamples) %>% 
    mutate(Braak_LBP = ifelse(asyn == 'yes', 6, 0)) %>% 
    dplyr::select(sample, age, sex, Braak_tau = braakbraak, APOE_WGS, Braak_LBP) %>% 
    left_join(., stab1_ad, by = 'sample') %>% 
    left_join(., stab1_dlb, by = 'sample') %>% 
    left_join(., stab1_pd, by = 'sample') %>% 
    left_join(., stab1_mapt, by = 'sample') %>% 
    left_join(., stab1_gba, by = 'sample') %>% 
    left_join(., stab1_snca, by = 'sample')
write_tsv(rzids_export, 'paper_export/STable_1.tsv')
```

# 4 Session Info

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
    ##  [1] GenomicFeatures_1.58.0      AnnotationDbi_1.68.0       
    ##  [3] data.table_1.16.4           rtracklayer_1.66.0         
    ##  [5] magrittr_2.0.3              VariantAnnotation_1.52.0   
    ##  [7] Rsamtools_2.22.0            Biostrings_2.74.1          
    ##  [9] XVector_0.46.0              SummarizedExperiment_1.36.0
    ## [11] Biobase_2.66.0              GenomicRanges_1.58.0       
    ## [13] GenomeInfoDb_1.42.3         IRanges_2.40.1             
    ## [15] S4Vectors_0.44.0            MatrixGenerics_1.18.1      
    ## [17] matrixStats_1.5.0           BiocGenerics_0.52.0        
    ## [19] patchwork_1.3.0             scales_1.3.0               
    ## [21] ggsci_3.2.0                 ggpubr_0.6.0               
    ## [23] broom_1.0.7                 lubridate_1.9.4            
    ## [25] forcats_1.0.0               stringr_1.5.1              
    ## [27] dplyr_1.1.4                 purrr_1.0.4                
    ## [29] readr_2.1.5                 tidyr_1.3.1                
    ## [31] tibble_3.2.1                ggplot2_3.5.1              
    ## [33] tidyverse_2.0.0            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3                bitops_1.0-9             rlang_1.1.5             
    ##  [4] compiler_4.4.2           RSQLite_2.3.9            png_0.1-8               
    ##  [7] vctrs_0.6.5              pkgconfig_2.0.3          crayon_1.5.3            
    ## [10] fastmap_1.2.0            backports_1.5.0          labeling_0.4.3          
    ## [13] utf8_1.2.4               rmarkdown_2.29           tzdb_0.4.0              
    ## [16] UCSC.utils_1.2.0         bit_4.5.0.1              xfun_0.50               
    ## [19] zlibbioc_1.52.0          cachem_1.1.0             jsonlite_1.8.9          
    ## [22] blob_1.2.4               DelayedArray_0.32.0      BiocParallel_1.38.0     
    ## [25] parallel_4.4.2           R6_2.6.0                 stringi_1.8.4           
    ## [28] car_3.1-3                knitr_1.49               Matrix_1.7-2            
    ## [31] timechange_0.3.0         tidyselect_1.2.1         rstudioapi_0.17.1       
    ## [34] abind_1.4-8              yaml_2.3.10              codetools_0.2-20        
    ## [37] curl_6.2.0               lattice_0.22-6           withr_3.0.2             
    ## [40] KEGGREST_1.46.0          evaluate_1.0.3           pillar_1.10.1           
    ## [43] carData_3.0-5            generics_0.1.3           vroom_1.6.5             
    ## [46] RCurl_1.98-1.16          hms_1.1.3                munsell_0.5.1           
    ## [49] glue_1.8.0               tools_4.4.2              BiocIO_1.16.0           
    ## [52] BSgenome_1.74.0          ggsignif_0.6.4           GenomicAlignments_1.42.0
    ## [55] XML_3.99-0.18            grid_4.4.2               colorspace_2.1-1        
    ## [58] GenomeInfoDbData_1.2.13  restfulr_0.0.15          Formula_1.2-5           
    ## [61] cli_3.6.4                S4Arrays_1.6.0           gtable_0.3.6            
    ## [64] rstatix_0.7.2            digest_0.6.37            SparseArray_1.6.1       
    ## [67] rjson_0.2.23             farver_2.1.2             memoise_2.0.1           
    ## [70] htmltools_0.5.8.1        lifecycle_1.0.4          httr_1.4.7              
    ## [73] bit64_4.6.0-1
