# Data and code for 'Sensory traps drive widespread bat-window collisions.'



# This repository contains:
## R/

Scripts used to read, tidy, standardize, and analyze data. Scripts are organized
into two tiers:

**Tier 1 — Reproducible from shared data:** Scripts 2–4 read only from
`data/derived/` (plus public APIs for map data in script 4) and reproduce all
Kansas City analyses and figures.

**Tier 2 — Data assembly pipeline:** Scripts 5–6 download and manually curate
iNaturalist records. They require raw data files, API access, and manual review
steps not included in this repository. `data/derived/iNaturalist records.csv`
is the shareable output of this pipeline.

Script 1 tidies the raw KC survey spreadsheets into `data/derived/`; it can be
skipped if starting from the shared derived files.

| Script | Description |
| ------ | ----------- |
| `0_funs.R` | Setup and functions used in subsequent scripts. |
| `1_setup and load KC data.R` | Load and tidy data from semi-structured surveys in Kansas City, MO ("KC"). The primary outputs are `data/derived/structured_surveys_schedule.csv` and `data/derived/structured_surveys_bats_discovered.csv`. |
| `2_plot and summarize KC records.R` | Plots and summaries from KC semi-structured surveys. |
| `3_model KC record timing.R` | Models predicting bat discoveries by day of year from KC surveys. |
| `4_map KC records.R` | Maps of KC survey records. Requires local NLCD data for the regional land cover map; see comments in script for directory structure. |
| `5_get iNat records.R` | Download candidate records from iNaturalist. After manual cleaning and removal of unlicensed data, the shareable output is `data/derived/iNaturalist records.csv`. |
| `6_analyze iNat records.R` | Analysis and figures from iNaturalist records. Requires `data/iNat_observations_tidy_manualChecks.csv` (not publicly archived due to licensing constraints). |


## data/derived/

Tidied data created or synthesized in the course of this study. **Scripts 2–4
can be run in full starting from the files in this directory.**

| File | Description |
| ---- | ----------- |
| `structured_surveys_schedule.csv` | Dates of semi-structured surveys in KC. Contains columns: <br> `survey` = Was survey conducted on specified date? (boolean) <br> `date` = Date of survey in yyyy-mm-dd format.<br> `yday` = Day of year of survey (1-365).<br>`yday_bin7` = Day of year of survey, divided into 7-day bins. |
| `structured_surveys_bats_discovered.csv` | Discovery of bats in KC. Contains columns: <br> `id` = ID of bat discovered (numeric)<br> `date` = Date of discovery in yyyy-mm-dd format.<br> `yday` = Day of year of discovery (1-365).<br> `species` = Species of discovered bat (common name)<br> `plotGroup` = Species group of discovered bat (most common species + Other category)<br> `locality` = Notes on locality of discovery<br> `Status` = Notes on status (typically alive vs. dead)<br> `Description Where Found` = Notes on specific conditions of discovery (e.g., on sidewalk)<br> `Notes` = Additional notes<br> `Building_side` = Cardinal direction of discovery relative to nearest building<br> `paired` = Was another bat found in immediate vicinity (Y/N). |
| `iNaturalist records.csv` | Retained iNaturalist records from data synthesis. Observations associated with restrictive licenses (e.g., "All rights reserved") are represented by ID and URL only. Contains columns: <br> `id` = iNaturalist ID number of observation<br> `url` = URL of observation <br> `notes` = Notes derived from manual checks of all observations. <br> `scientific_name` = Scientific name of organism from observation. <br> `common_name` = Common name of organism from observation<br> `description` = Text description provided by observer of observation. <br>  `user_login` = Username of user who made observation <br> `observed_on` = Datetime of observation <br> `license` = License of observation (at time of download)<br> `family_name` = Scientific name of family <br> `adm0_a3` = Three-letter country code where record was made <br>  `continent` = Continent of record |


## tmp/

Intermediate plot objects cached between scripts (gitignored; created on first
run). Scripts 2 and 3 use this directory to pass assembled ggplot objects between
steps when building multi-panel figures.


## out/

Core outputs of analyses.

| Directory | Description |
| --------- | ----------- |
| `figs/` | Figures used in manuscript and SI. |
| `models/` | Fitted model objects generated in `3_model KC record timing.R` (`.rds` format, cached by brms). |


## Session info
Code was run most recently with the following session information:
R version 4.4.3 (2025-02-28)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.8.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets 
[7] methods   base     

other attached packages:
 [1] tidybayes_3.0.7     rstan_2.32.7        StanHeaders_2.32.10
 [4] brms_2.23.0         Rcpp_1.1.0          sjPlot_2.9.0       
 [7] glmmTMB_1.1.13      magick_2.9.0        rinat_0.1.10       
[10] geodata_0.6-6       data.table_1.17.8   patchwork_1.3.2    
[13] rnaturalearth_1.1.0 terra_1.8-70        sf_1.0-23          
[16] tidyterra_0.7.2     readxl_1.4.5        lubridate_1.9.4    
[19] forcats_1.0.1       stringr_1.5.2       dplyr_1.1.4        
[22] purrr_1.1.0         readr_2.1.5         tidyr_1.3.1        
[25] tibble_3.3.0        ggplot2_4.0.0       tidyverse_2.0.0    

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      tensorA_0.36.2.1       
  [3] rstudioapi_0.17.1       jsonlite_2.0.0         
  [5] datawizard_1.3.0        magrittr_2.0.4         
  [7] estimability_1.5.1      farver_2.1.2           
  [9] nloptr_2.2.1            ragg_1.5.0             
 [11] vctrs_0.6.5             minqa_1.2.8            
 [13] distributional_0.5.0    curl_7.0.0             
 [15] cellranger_1.1.0        sjmisc_2.8.11          
 [17] KernSmooth_2.23-26      plyr_1.8.9             
 [19] sandwich_3.1-1          emmeans_1.11.2-8       
 [21] zoo_1.8-14              TMB_1.9.18             
 [23] commonmark_2.0.0        lifecycle_1.0.4        
 [25] pkgconfig_2.0.3         sjlabelled_1.2.0       
 [27] Matrix_1.7-4            R6_2.6.1               
 [29] rbibutils_2.3           numDeriv_2016.8-1.1    
 [31] colorspace_2.1-2        ps_1.9.1               
 [33] textshaping_1.0.4       labeling_0.4.3         
 [35] timechange_0.3.0        httr_1.4.7             
 [37] abind_1.4-8             mgcv_1.9-3             
 [39] compiler_4.4.3          proxy_0.4-27           
 [41] bit64_4.6.0-1           withr_3.0.2            
 [43] inline_0.3.21           S7_0.2.0               
 [45] backports_1.5.0         DBI_1.2.3              
 [47] QuickJSR_1.8.1          pkgbuild_1.4.8         
 [49] maps_3.4.3              MASS_7.3-65            
 [51] rappdirs_0.3.3          classInt_0.4-11        
 [53] loo_2.8.0               tools_4.4.3            
 [55] units_1.0-0             glue_1.8.0             
 [57] rnaturalearthdata_1.0.0 nlme_3.1-168           
 [59] gridtext_0.1.5          grid_4.4.3             
 [61] cmdstanr_0.8.1          checkmate_2.3.3        
 [63] reshape2_1.4.4          generics_0.1.4         
 [65] gtable_0.3.6            tzdb_0.5.0             
 [67] class_7.3-23            hms_1.1.4              
 [69] xml2_1.4.0              utf8_1.2.6             
 [71] ggdist_3.3.3            pillar_1.11.1          
 [73] markdown_2.0            vroom_1.6.6            
 [75] posterior_1.6.1         ggtext_0.1.2           
 [77] lattice_0.22-7          bit_4.6.0              
 [79] tidyselect_1.2.1        knitr_1.50             
 [81] arrayhelpers_1.1-0      gridExtra_2.3          
 [83] reformulas_0.4.1        V8_8.0.1               
 [85] litedown_0.7            svglite_2.2.2          
 [87] stats4_4.4.3            xfun_0.53              
 [89] bridgesampling_1.1-2    matrixStats_1.5.0      
 [91] stringi_1.8.7           rematch_2.0.0          
 [93] boot_1.3-32             evaluate_1.0.5         
 [95] codetools_0.2-20        cli_3.6.5              
 [97] RcppParallel_5.1.11-1   xtable_1.8-4           
 [99] systemfonts_1.3.1       Rdpack_2.6.4           
[101] processx_3.8.6          dichromat_2.0-0.1      
[103] MCMCvis_0.16.3          coda_0.19-4.1          
[105] svUnit_1.0.8            parallel_4.4.3         
[107] rstantools_2.5.0        bayestestR_0.17.0      
[109] bayesplot_1.14.0        Brobdingnag_1.2-9      
[111] lme4_1.1-37             viridisLite_0.4.2      
[113] mvtnorm_1.3-3           scales_1.4.0           
[115] e1071_1.7-16            crayon_1.5.3           
[117] insight_1.4.2           rlang_1.1.6         
