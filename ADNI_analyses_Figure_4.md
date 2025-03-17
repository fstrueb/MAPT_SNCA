ADNI tau PET analysis
================
<felix.struebing@med.uni-muenchen.de>
2025-03-17

- [0.1 Load ADNI data](#01-load-adni-data)
- [0.2 Effect of SAA on tau PET
  (SUVR)](#02-effect-of-saa-on-tau-pet-suvr)
- [0.3 Session info](#03-session-info)

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggpubr)
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
source('ggplot_theme_patchwork.R')
```

## 0.1 Load ADNI data

``` r
mytable <- read_csv('resources/ADNI_DF_11_2022_22_11_2024.csv') %>% 
    mutate(saa = CSF_AlphaSyn_seeding)
```

    ## New names:
    ## Rows: 9491 Columns: 1675
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (49): ID, Phase, DX, MCI_memory_features, DX_Dementia_severity, DX_De... dbl
    ## (1586): ...1, RID, SITEID, year_of_cognitive_symptom_onset, education_y... lgl
    ## (12): amyloid.SUVR.DK.lcorpuscallosum, amyloid.SUVR.DK.rcorpuscallosu... date
    ## (28): EXAMDATE_diagnosis, birthdate_full, USERDATE, EXAMDATE_amyloid,...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
dim(mytable)
```

    ## [1] 9491 1676

## 0.2 Effect of SAA on tau PET (SUVR)

``` r
model <- lm(tau.SUVR.DK.lsuperiorfrontal ~ saa + sex + age, data = mytable)
summary(model)
```

    ## 
    ## Call:
    ## lm(formula = tau.SUVR.DK.lsuperiorfrontal ~ saa + sex + age, 
    ##     data = mytable)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.32569 -0.10973 -0.04196  0.03715  1.17410 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.290080   0.102759  12.554  < 2e-16 ***
    ## saapositive  0.096726   0.029292   3.302  0.00105 ** 
    ## sexmale     -0.016123   0.022119  -0.729  0.46651    
    ## age         -0.002902   0.001428  -2.032  0.04287 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2061 on 365 degrees of freedom
    ##   (9122 observations deleted due to missingness)
    ## Multiple R-squared:  0.03585,    Adjusted R-squared:  0.02792 
    ## F-statistic: 4.524 on 3 and 365 DF,  p-value: 0.00395

``` r
# Get the estimated marginal means (adjusted for sex and age)
emm_results <- emmeans(model, ~ saa, type = "response")
emm_df <- as.data.frame(emm_results)
sfig_ADNI <- ggplot(emm_df, aes(x = saa, y = emmean, fill = saa)) +
    geom_col(color = 'black') +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    geom_signif(comparisons = list(c('negative', 'positive')), annotations = '***', y_position = 1.35, tip_length = 0.1) +
    labs(x = "aSyn SAA status", y = "Adjusted tau.SUVR (SFG)") +
    labs(title = 'SFG tau PET by aSyn SAA') +
    guides(fill = 'none') +
    ylim(0, 1.4) +
    coord_flip() +
    theme_patch()
```

``` r
sfig_ADNI
```

![](gh_v2_ADNI_genotyping_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# pdf('SFig_ADNI.pdf', width = 9, height = 4)
# sfig_ADNI
# dev.off()
```

Fit nested linear models (adjusted for sex and age) for all tau SUVR
Schaefer regions and correct the resulting p-value with
Benjamini-Hochberg’s procedure:

``` r
adni_nested <- mytable %>% 
    dplyr::select(sex, saa, age, starts_with('tau.SUVR.Schaefer200.ROI.idx')) %>% 
    filter(!is.na(saa) & !is.na(tau.SUVR.Schaefer200.ROI.idx.1)) %>% 
    pivot_longer(-c(sex, age, saa)) %>% 
    group_by(name) %>% 
    nest() %>% 
    mutate(fit = map(data, ~lm(value ~ saa + sex + age, data = .)),
           tidy = map(fit, broom::tidy)) %>% 
    unnest(tidy) %>% 
    dplyr::select(name, term, estimate, std.error, p.value) %>% 
    ungroup() %>% 
    filter(term == 'saapositive') %>% 
    mutate(padj = p.adjust(p.value, method = 'BH')) %>% 
    arrange(padj)
```

Number of significant results:

``` r
table(adni_nested$padj < 0.05)
```

    ## 
    ## FALSE  TRUE 
    ##     6   194

Export as supplementary figure:

``` r
adni_supp <- adni_nested %>% 
    mutate(significant = ifelse(padj < 0.05, 'yes', 'no')) %>% 
    ggplot(aes(x = reorder(name, estimate), y = estimate, ymax = estimate + std.error, ymin = estimate - std.error, color = significant)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_errorbar(width = 0.2) +
    labs(title = 'Tau PET signal: Effect of aSyn SAA', y = 'Estimate (Schaefer ROI SUVR ~ SAA + age + sex)', x = '', color = 'FDR < 0.05') +
    coord_flip() +
    theme_patch()
# adni_supp
```

``` r
# pdf('SFig_ADNI_allres.pdf', height = 26, width = 10)
# adni_supp
# dev.off()
```

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] emmeans_1.10.7  ggpubr_0.6.0    lubridate_1.9.4 forcats_1.0.0  
    ##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5    
    ##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       xfun_0.50          rstatix_0.7.2      lattice_0.22-6    
    ##  [5] tzdb_0.4.0         vctrs_0.6.5        tools_4.4.2        generics_0.1.3    
    ##  [9] parallel_4.4.2     sandwich_3.1-1     pkgconfig_2.0.3    Matrix_1.7-2      
    ## [13] lifecycle_1.0.4    farver_2.1.2       compiler_4.4.2     munsell_0.5.1     
    ## [17] codetools_0.2-20   carData_3.0-5      htmltools_0.5.8.1  yaml_2.3.10       
    ## [21] Formula_1.2-5      crayon_1.5.3       pillar_1.10.1      car_3.1-3         
    ## [25] MASS_7.3-64        abind_1.4-8        multcomp_1.4-28    tidyselect_1.2.1  
    ## [29] digest_0.6.37      mvtnorm_1.3-3      stringi_1.8.4      labeling_0.4.3    
    ## [33] splines_4.4.2      fastmap_1.2.0      grid_4.4.2         colorspace_2.1-1  
    ## [37] cli_3.6.4          magrittr_2.0.3     survival_3.8-3     broom_1.0.7       
    ## [41] TH.data_1.1-3      withr_3.0.2        scales_1.3.0       backports_1.5.0   
    ## [45] bit64_4.6.0-1      estimability_1.5.1 timechange_0.3.0   rmarkdown_2.29    
    ## [49] bit_4.5.0.1        ggsignif_0.6.4     zoo_1.8-12         hms_1.1.3         
    ## [53] coda_0.19-4.1      evaluate_1.0.3     knitr_1.49         rlang_1.1.5       
    ## [57] xtable_1.8-4       glue_1.8.0         rstudioapi_0.17.1  vroom_1.6.5       
    ## [61] R6_2.6.0
