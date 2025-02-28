
# Systems analysis unravels a common rural-urban gradient in immunological profile, function and metabolic dependencies

```
Author      : Mikhael D. Manurung
Email       : m.d.manurung@lumc.nl
Affiliation : Leiden University Center for Infectious Diseases (LU-CID)
```

This repository contains data and analysis scripts to reproduce figures shown in the manuscript.

Folders:

```
data/       # processed data for downstream analyses
outputs/    # models, plots, and tables are saved here
src/        # analysis scripts are saved here
```

Note that the following data objects have to be downloaded from Zenodo:

```
# Zenodo URL: https://zenodo.org/records/14852065

data/met_inhibitions/integrated_allLineages.qs
data/met_inhibitions/FigS8_B_umap_coordinates.qs
data/EXV_sce_annotated.qs
data/PMA_data.qs
data/MPL_data.qs
```

Purpose of scripts in `src/` folder:

```
Fig1_BC.R
# script to plot age and sex distribution as well as IgG1 galactosylation
# outputs: Fig1_B, Fig1_C


Fig2_train_mofa.R       
# script to run MOFA and select final MOFA model
# outputs: Fig1_B


Fig2_principal_curve.R  
# script to plot MDS and perform principal curve analysis
# depends on: Fig2_train_mofa.R
# outputs: Fig1_CD, FigS2_A


FigS1_ABC.R
# script to generate figures for Fig S1A, B, and C
# outputs: FigS1_A-C


FigS4_compare_CAA_rural.R
# outputs: FigS4 and TableS6

Fig8_B.R
# outputs: Fig8_B and FigS8_C

Fig8_C.R
# outputs: Fig8_C

Fig8_D.R
# outputs: Fig8_D and FigS8_D

FigS8_AB_calculate_iMFI.R
# umap and marker heatmaps of metabolic enzyme inhibition spectral flow data
# calculate iMFI for use in Fig8
# outputs: FigS8_AB

```



`sessionInfo()` outputs:

```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Dutch_Netherlands.utf8  LC_CTYPE=Dutch_Netherlands.utf8    LC_MONETARY=Dutch_Netherlands.utf8 LC_NUMERIC=C                      
[5] LC_TIME=Dutch_Netherlands.utf8    

time zone: Europe/Amsterdam
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggnewscale_0.5.0 vegan_2.6-4      lattice_0.21-8   permute_0.9-7    ggbeeswarm_0.7.2 ggpubr_0.6.0     princurve_2.1.6  PMCMRplus_1.9.10
 [9] emmeans_1.10.4   nlme_3.1-162     Hmisc_5.1-3      patchwork_1.3.0  glue_1.7.0       here_1.0.1       lubridate_1.9.3  forcats_1.0.0   
[17] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.2      readr_2.1.5      tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.1    tidyverse_2.0.0 
[25] MOFA2_1.12.1    

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3    rstudioapi_0.16.0     jsonlite_1.8.8        tidydr_0.0.5          magrittr_2.0.3        TH.data_1.1-2        
  [7] estimability_1.5.1    SuppDists_1.1-9.7     farver_2.1.1          corrplot_0.94         rmarkdown_2.28        fs_1.6.4             
 [13] ragg_1.3.3            zlibbioc_1.48.2       vctrs_0.6.5           memoise_2.0.1         base64enc_0.1-3       rstatix_0.7.2        
 [19] htmltools_0.5.8.1     S4Arrays_1.2.1        broom_1.0.7           Rhdf5lib_1.24.2       BWStest_0.2.3         SparseArray_1.2.4    
 [25] Formula_1.2-5         rhdf5_2.44.0          htmlwidgets_1.6.4     basilisk_1.14.3       plyr_1.8.9            sandwich_3.1-1       
 [31] zoo_1.8-12            cachem_1.0.8          conflicted_1.2.0      lifecycle_1.0.4       pkgconfig_2.0.3       Matrix_1.6-5         
 [37] R6_2.5.1              fastmap_1.1.1         MatrixGenerics_1.14.0 digest_0.6.35         colorspace_2.1-0      S4Vectors_0.40.2     
 [43] rprojroot_2.0.4       textshaping_0.4.0     filelock_1.0.3        labeling_0.4.3        fansi_1.0.6           timechange_0.3.0     
 [49] mgcv_1.8-42           abind_1.4-8           compiler_4.3.1        withr_3.0.1           htmlTable_2.4.3       backports_1.4.1      
 [55] carData_3.0-5         HDF5Array_1.28.1      ggsignif_0.6.4        MASS_7.3-60           DelayedArray_0.28.0   tools_4.3.1          
 [61] vipor_0.4.7           foreign_0.8-84        beeswarm_0.4.0        nnet_7.3-19           rhdf5filters_1.12.1   grid_4.3.1           
 [67] checkmate_2.3.2       Rtsne_0.17            cluster_2.1.4         reshape2_1.4.4        generics_0.1.3        gtable_0.3.5         
 [73] tzdb_0.4.0            data.table_1.15.4     hms_1.1.3             car_3.1-2             utf8_1.2.4            XVector_0.42.0       
 [79] BiocGenerics_0.48.1   ggrepel_0.9.5         pillar_1.9.0          yulab.utils_0.1.7     splines_4.3.1         survival_3.5-5       
 [85] gmp_0.7-5             tidyselect_1.2.1      knitr_1.48            gridExtra_2.3         IRanges_2.36.0        stats4_4.3.1         
 [91] xfun_0.47             matrixStats_1.2.0     pheatmap_1.0.12       stringi_1.8.3         ggfun_0.1.6           evaluate_1.0.0       
 [97] codetools_0.2-19      kSamples_1.2-10       multcompView_0.1-10   cli_3.6.2             uwot_0.2.2            rpart_4.1.19         
[103] systemfonts_1.1.0     xtable_1.8-4          reticulate_1.36.1     munsell_0.5.1         Rcpp_1.0.12           dir.expiry_1.10.0    
[109] png_0.1-8             parallel_4.3.1        basilisk.utils_1.14.1 Rmpfr_0.9-5           viridisLite_0.4.2     mvtnorm_1.2-4        
[115] scales_1.3.0          crayon_1.5.3          rlang_1.1.3           cowplot_1.1.3         multcomp_1.4-26     
```
