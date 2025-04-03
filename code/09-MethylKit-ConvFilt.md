09-MethylKit, filtered for \>90% conversion efficiency
================
Zoe Dellaert
2025-02-09

- [0.1 MethylKit - Reads filtered for \>90% conversion
  efficiency](#01-methylkit---reads-filtered-for-90-conversion-efficiency)
- [0.2 Managing Packages Using Renv](#02-managing-packages-using-renv)
- [0.3 Load packages](#03-load-packages)
  - [0.3.1 batch effects](#031-batch-effects)
  - [0.3.2 other possible filtering](#032-other-possible-filtering)
  - [0.3.3 Identify DML](#033-identify-dml)
- [0.4 Further look at genome wide
  methylation](#04-further-look-at-genome-wide-methylation)
  - [0.4.1 Annotation](#041-annotation)
- [0.5 Are any DMGs DMLs?](#05-are-any-dmgs-dmls)
  - [0.5.1 Read counts](#051-read-counts)

## 0.1 MethylKit - Reads filtered for \>90% conversion efficiency

I am identifying differentially methylated loci using methylkit based on
[Yaamini Venkataraman’s code](https://osf.io/u46xj)

## 0.2 Managing Packages Using Renv

To run this code in my project using the renv environment, run the
following lines of code

``` r
install.packages("renv") #install the package on the new computer (may not be necessary if renv bootstraps itself as expected)
renv::restore() #reinstall all the package versions in the renv lockfile
```

## 0.3 Load packages

``` r
require("methylKit")
```

    ## Loading required package: methylKit

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

``` r
require("tidyverse")
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks BiocGenerics::combine()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::select()       masks methylKit::select()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ✖ tidyr::unite()        masks methylKit::unite()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
require("vegan")
```

    ## Loading required package: vegan
    ## Loading required package: permute
    ## Loading required package: lattice

``` r
require("gplots")
```

    ## Loading required package: gplots
    ## 
    ## Attaching package: 'gplots'
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     space
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     space
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
require("ggplot2")
require("ggrepel")
```

    ## Loading required package: ggrepel

``` r
require("ggpmisc")
```

    ## Loading required package: ggpmisc
    ## Loading required package: ggpp
    ## Registered S3 methods overwritten by 'ggpp':
    ##   method                  from   
    ##   heightDetails.titleGrob ggplot2
    ##   widthDetails.titleGrob  ggplot2
    ## 
    ## Attaching package: 'ggpp'
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     annotate

``` r
require("dichromat")
```

    ## Loading required package: dichromat

``` r
require("readr")
require("genomationData")
```

    ## Loading required package: genomationData

``` r
require("genomation")
```

    ## Loading required package: genomation
    ## Loading required package: grid

    ## Warning: replacing previous import 'Biostrings::pattern' by 'grid::pattern'
    ## when loading 'genomation'

    ## 
    ## Attaching package: 'genomation'
    ## 
    ## The following objects are masked from 'package:methylKit':
    ## 
    ##     getFeatsWithTargetsStats, getFlanks, getMembers,
    ##     getTargetAnnotationStats, plotTargetAnnotation

``` r
sessionInfo() #provides list of loaded packages and version of R.
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices datasets  utils    
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] genomation_1.38.0     genomationData_1.38.0 dichromat_2.0-0.1    
    ##  [4] ggpmisc_0.6.1         ggpp_0.5.8-1          ggrepel_0.9.6        
    ##  [7] gplots_3.2.0          vegan_2.6-10          lattice_0.22-6       
    ## [10] permute_0.9-7         lubridate_1.9.4       forcats_1.0.0        
    ## [13] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.4          
    ## [16] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
    ## [19] ggplot2_3.5.1         tidyverse_2.0.0       methylKit_1.32.0     
    ## [22] GenomicRanges_1.58.0  GenomeInfoDb_1.42.3   IRanges_2.40.1       
    ## [25] S4Vectors_0.44.0      BiocGenerics_0.52.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.17.1           jsonlite_1.8.9             
    ##  [3] magrittr_2.0.3              rmarkdown_2.29             
    ##  [5] BiocIO_1.16.0               zlibbioc_1.52.0            
    ##  [7] vctrs_0.6.5                 Rsamtools_2.22.0           
    ##  [9] RCurl_1.98-1.16             htmltools_0.5.8.1          
    ## [11] S4Arrays_1.6.0              polynom_1.4-1              
    ## [13] plotrix_3.8-4               curl_6.2.0                 
    ## [15] SparseArray_1.6.1           KernSmooth_2.23-26         
    ## [17] plyr_1.8.9                  impute_1.80.0              
    ## [19] GenomicAlignments_1.42.0    lifecycle_1.0.4            
    ## [21] pkgconfig_2.0.3             Matrix_1.7-2               
    ## [23] R6_2.6.0                    fastmap_1.2.0              
    ## [25] GenomeInfoDbData_1.2.13     MatrixGenerics_1.18.1      
    ## [27] digest_0.6.37               numDeriv_2016.8-1.1        
    ## [29] colorspace_2.1-1            timechange_0.3.0           
    ## [31] httr_1.4.7                  abind_1.4-8                
    ## [33] mgcv_1.9-1                  compiler_4.4.0             
    ## [35] withr_3.0.2                 BiocParallel_1.40.0        
    ## [37] R.utils_2.12.3              MASS_7.3-64                
    ## [39] quantreg_6.00               DelayedArray_0.32.0        
    ## [41] rjson_0.2.23                gtools_3.9.5               
    ## [43] caTools_1.18.3              tools_4.4.0                
    ## [45] R.oo_1.27.0                 glue_1.8.0                 
    ## [47] restfulr_0.0.15             nlme_3.1-167               
    ## [49] gridBase_0.4-7              cluster_2.1.8              
    ## [51] reshape2_1.4.4              generics_0.1.3             
    ## [53] gtable_0.3.6                BSgenome_1.74.0            
    ## [55] tzdb_0.4.0                  R.methodsS3_1.8.2          
    ## [57] seqPattern_1.38.0           data.table_1.16.4          
    ## [59] hms_1.1.3                   XVector_0.46.0             
    ## [61] pillar_1.10.1               fastseg_1.52.0             
    ## [63] emdbook_1.3.13              limma_3.62.2               
    ## [65] splines_4.4.0               renv_1.1.1                 
    ## [67] survival_3.8-3              rtracklayer_1.66.0         
    ## [69] SparseM_1.84-2              tidyselect_1.2.1           
    ## [71] Biostrings_2.74.1           knitr_1.49                 
    ## [73] SummarizedExperiment_1.36.0 xfun_0.50                  
    ## [75] Biobase_2.66.0              statmod_1.5.0              
    ## [77] matrixStats_1.5.0           stringi_1.8.4              
    ## [79] UCSC.utils_1.2.0            yaml_2.3.10                
    ## [81] evaluate_1.0.3              codetools_0.2-20           
    ## [83] bbmle_1.0.25.1              qvalue_2.38.0              
    ## [85] BiocManager_1.30.25         cli_3.6.4                  
    ## [87] munsell_0.5.1               Rcpp_1.0.14                
    ## [89] coda_0.19-4.1               bdsmatrix_1.3-7            
    ## [91] XML_3.99-0.18               parallel_4.4.0             
    ## [93] MatrixModels_0.5-3          mclust_6.1.1               
    ## [95] bitops_1.0-9                mvtnorm_1.3-3              
    ## [97] scales_1.3.0                crayon_1.5.3               
    ## [99] rlang_1.1.5

``` r
meta <- read.csv("../data_WGBS/LCM_WGBS_metadata.csv", sep = ",", header = TRUE) %>%
  mutate(Section_Date = as.character(Section_Date), LCM_Date = as.character(LCM_Date),DNA_Extraction_Date = as.character(DNA_Extraction_Date))

meta <- meta %>% arrange(Sample)
```

``` r
#file_list <- list.files("/scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/deduplicated/min_efficiency_test",pattern = "^min_90.*_CpG.methylKit$",  full.names = TRUE, include.dirs = FALSE)

file_list <- list.files("../output_WGBS/methylseq_V3_bwa_test/methyldackel/min_efficiency_test_new",pattern = "^min_90.*_CpG.methylKit$",  full.names = TRUE, include.dirs = FALSE)

sample <- gsub("_CpG.methylKit", "", basename(file_list) )
sample <- gsub("min_90_", "", sample)

sample == meta$Sample #the files and metadata are in the same order
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
tissue <- meta$Tissue
tissue_binary <- gsub("Aboral", "1", tissue)
tissue_binary <- gsub("OralEpi", "0", tissue_binary)
tissue_binary <- as.numeric(tissue_binary)
fragment <- meta$Fragment

methylObj=methRead(as.list(file_list),
           sample.id = as.list(sample),
           assembly = "Pacuta",
           treatment = tissue_binary,
           context = "CpG",
           mincov = 5
           )
```

    ## Received list of locations.

    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.
    ## Reading file.

``` r
# save(methylObj, file = "../output_WGBS/MethylKit.RData")
```

``` r
#load("../output_WGBS/MethylKit.RData")

getMethylationStats(methylObj[[2]],plot=FALSE,both.strands=FALSE)
```

    ## methylation statistics per base
    ## summary:
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   0.000   0.000   0.000   2.573   0.000 100.000 
    ## percentiles:
    ##    0%   10%   20%   30%   40%   50%   60%   70%   80%   90%   95%   99% 99.5% 
    ##     0     0     0     0     0     0     0     0     0     0     0   100   100 
    ## 99.9%  100% 
    ##   100   100

``` r
getMethylationStats(methylObj[[2]],plot=TRUE,both.strands=FALSE)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
getCoverageStats(methylObj[[2]],plot=TRUE,both.strands=FALSE)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
getCoverageStats(methylObj[[2]],plot=TRUE,both.strands=TRUE)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
filtered_methylObj=filterByCoverage(methylObj,lo.count=5,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=99.9)

filtered_methylObj_norm <- filtered_methylObj %>% methylKit::normalizeCoverage(.)
```

``` r
meth_filter=methylKit::unite(filtered_methylObj_norm)
```

    ## uniting...

``` r
meth_filter_destrand=methylKit::unite(filtered_methylObj_norm,destrand = TRUE)
```

    ## destranding...
    ## uniting...

``` r
clusterSamples(meth_filter, dist="correlation", method="ward", plot=TRUE)
```

    ## The "ward" method has been renamed to "ward.D"; note new "ward.D2"

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

    ## 
    ## Call:
    ## hclust(d = d, method = HCLUST.METHODS[hclust.method])
    ## 
    ## Cluster method   : ward.D 
    ## Distance         : pearson 
    ## Number of objects: 10

``` r
clusterSamples(meth_filter_destrand, dist="correlation", method="ward", plot=TRUE)
```

    ## The "ward" method has been renamed to "ward.D"; note new "ward.D2"

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

    ## 
    ## Call:
    ## hclust(d = d, method = HCLUST.METHODS[hclust.method])
    ## 
    ## Cluster method   : ward.D 
    ## Distance         : pearson 
    ## Number of objects: 10

``` r
PCASamples(meth_filter)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
PCASamples(meth_filter_destrand)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
getCorrelation(meth_filter_destrand,plot=TRUE)
```

    ##            LCM_1    LCM_11    LCM_12    LCM_17    LCM_18    LCM_24    LCM_25
    ## LCM_1  1.0000000 0.6548849 0.6442972 0.6675210 0.7071313 0.6183485 0.6590492
    ## LCM_11 0.6548849 1.0000000 0.6131597 0.6221372 0.6138284 0.5603577 0.6079632
    ## LCM_12 0.6442972 0.6131597 1.0000000 0.6808556 0.6822028 0.6467358 0.6924756
    ## LCM_17 0.6675210 0.6221372 0.6808556 1.0000000 0.6908618 0.6052214 0.6836569
    ## LCM_18 0.7071313 0.6138284 0.6822028 0.6908618 1.0000000 0.6289884 0.6983488
    ## LCM_24 0.6183485 0.5603577 0.6467358 0.6052214 0.6289884 1.0000000 0.6871517
    ## LCM_25 0.6590492 0.6079632 0.6924756 0.6836569 0.6983488 0.6871517 1.0000000
    ## LCM_3  0.5934479 0.5942247 0.5873359 0.6078850 0.5753703 0.6234846 0.6210878
    ## LCM_32 0.6273917 0.5553906 0.5911450 0.5960248 0.5920193 0.5671652 0.6003708
    ## LCM_33 0.6219982 0.5252421 0.5934366 0.6355371 0.6420074 0.5373079 0.5931981
    ##            LCM_3    LCM_32    LCM_33
    ## LCM_1  0.5934479 0.6273917 0.6219982
    ## LCM_11 0.5942247 0.5553906 0.5252421
    ## LCM_12 0.5873359 0.5911450 0.5934366
    ## LCM_17 0.6078850 0.5960248 0.6355371
    ## LCM_18 0.5753703 0.5920193 0.6420074
    ## LCM_24 0.6234846 0.5671652 0.5373079
    ## LCM_25 0.6210878 0.6003708 0.5931981
    ## LCM_3  1.0000000 0.5827830 0.4866721
    ## LCM_32 0.5827830 1.0000000 0.5748508
    ## LCM_33 0.4866721 0.5748508 1.0000000

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### 0.3.1 batch effects

``` r
as=assocComp(mBase=meth_filter_destrand,dplyr::select(meta,c("PCR_ReAmp_Cycles", "Fragment")))
as
```

    ## $pcs
    ##               PC1         PC2         PC3         PC4        PC5         PC6
    ## LCM_1  -0.3278981  0.12625480  0.18982046  0.18376942  0.2390968 -0.38714766
    ## LCM_11 -0.3052206 -0.18122081  0.38587520  0.67073994  0.1138058 -0.09699804
    ## LCM_12 -0.3250024  0.01216415 -0.29540627  0.07384750  0.1930995  0.53779603
    ## LCM_17 -0.3278144  0.14561919 -0.06870506  0.13823951 -0.3381564  0.43750457
    ## LCM_18 -0.3300652  0.22475300 -0.18624015  0.14059548  0.1147250 -0.08109501
    ## LCM_24 -0.3118345 -0.33403137 -0.43656180 -0.27582151  0.1860107 -0.45649863
    ## LCM_25 -0.3305895 -0.08956165 -0.34202501 -0.05199291  0.1114517  0.09241968
    ## LCM_3  -0.3011445 -0.58621666  0.19818948 -0.14885317 -0.6119190 -0.02570500
    ## LCM_32 -0.3017975  0.01990903  0.58262496 -0.57650744  0.3888855  0.23218534
    ## LCM_33 -0.2983141  0.64475371  0.04582194 -0.19674136 -0.4386156 -0.28634537
    ##                PC7         PC8         PC9        PC10
    ## LCM_1   0.43704976  0.23049787 -0.30679726  0.51129504
    ## LCM_11 -0.43346114 -0.17856653  0.10903110 -0.14555043
    ## LCM_12 -0.30165952  0.59787664 -0.02718236  0.16062784
    ## LCM_17  0.19912895 -0.41379169 -0.56571458 -0.10408334
    ## LCM_18  0.46881302  0.16557232  0.33598091 -0.63739258
    ## LCM_24 -0.28432815 -0.09289277 -0.37400164 -0.22842675
    ## LCM_25  0.08723227 -0.51112209  0.52498405  0.44449806
    ## LCM_3   0.16824506  0.26897416  0.16555857  0.03250943
    ## LCM_32 -0.02318647 -0.11528313  0.01596961 -0.13222861
    ## LCM_33 -0.39162036  0.05812346  0.13108194  0.07867609
    ## 
    ## $vars
    ##  [1] 65.682286  5.671953  4.790223  4.629650  3.761410  3.578142  3.355369
    ##  [8]  2.984270  2.882112  2.664585
    ## 
    ## $association
    ##                        PC1        PC2       PC3       PC4       PC5       PC6
    ## PCR_ReAmp_Cycles 0.7968283 0.03306508 0.5372258 0.3863827 0.6475336 0.5721282
    ## Fragment         0.2854370 0.22495865 0.1991483 0.1991483 0.9689278 0.7424472
    ##                        PC7       PC8       PC9      PC10
    ## PCR_ReAmp_Cycles 0.6138154 0.5804624 0.8785944 0.1338345
    ## Fragment         0.1058444 0.4311711 0.9944651 0.6044507

### 0.3.2 other possible filtering

``` r
# get percent methylation matrix
pm=percMethylation(meth_filter_destrand)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# keep only CpG with standard deviations larger than 2%
meth <- meth_filter_destrand[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth_filter_destrand)
```

    ## [1] 2987

``` r
nrow(meth)
```

    ## [1] 1451

### 0.3.3 Identify DML

``` r
DMLStats_Tissue <- methylKit::calculateDiffMeth(meth_filter_destrand, overdispersion = "MN", test = "Chisq", mc.cores = 8) #Calculate differential methylation statistics and include covariate information.
```

    ## two groups detected:
    ##  will calculate methylation difference as the difference of
    ## treatment (group: 1) - control (group: 0)

``` r
head(DMLStats_Tissue) #Look at differential methylation output
```

    ##                                  chr   start     end strand     pvalue
    ## 1 Pocillopora_acuta_HIv2___Sc0000000  844614  844614      + 0.51009168
    ## 2 Pocillopora_acuta_HIv2___Sc0000000  844624  844624      + 0.97976035
    ## 3 Pocillopora_acuta_HIv2___Sc0000000  844628  844628      + 0.34743888
    ## 4 Pocillopora_acuta_HIv2___Sc0000000  867803  867803      + 0.22663944
    ## 5 Pocillopora_acuta_HIv2___Sc0000000 5859931 5859931      + 1.00000000
    ## 6 Pocillopora_acuta_HIv2___Sc0000000 8087338 8087338      + 0.01061017
    ##      qvalue   meth.diff
    ## 1 0.9453423  1.03812117
    ## 2 0.9453423  0.03189793
    ## 3 0.7517825 -3.71219844
    ## 4 0.6829760 -3.48837209
    ## 5 0.9453423  0.00000000
    ## 6 0.3947417 -8.13953488

``` r
# Filter DMRs with q-value < 0.05
significant_dmg <- getData(DMLStats_Tissue[DMLStats_Tissue$qvalue < 0.1, ])

# Create a data frame for plotting
plot_data <- data.frame(
  chr = significant_dmg$chr,
  start = significant_dmg$start,
  meth.diff = significant_dmg$meth.diff
)

# Count the number of positive and negative methylation differences
positive_count <- sum(significant_dmg$meth.diff > 0)
negative_count <- sum(significant_dmg$meth.diff < 0)

# Plot with counts added to the quadrants
ggplot(plot_data, aes(x = start, y = meth.diff)) +
  geom_point(alpha = 0.5) +  # Set alpha to reduce point transparency
  theme_minimal() +
  labs(title = "Significant Differentially Methylated Regions (q-value < 0.05)",
       x = "Genomic Position (start)",
       y = "Methylation Difference (%)") +
  theme(legend.position = "none") +  # Remove the legend
  # Add the count of positive and negative methylation differences as text annotations
  annotate("text", x = Inf, y = Inf, label = paste("Positive:", positive_count), 
           hjust = 1.1, vjust = 1.1, size = 4, color = "blue") +
  annotate("text", x = Inf, y = -Inf, label = paste("Negative:", negative_count), 
           hjust = 1.1, vjust = -0.1, size = 4, color = "red")
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
dml_df <- as.data.frame(DMLStats_Tissue)

# Volcano plot
ggplot(dml_df, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Differentially Methylated Loci (DMLs)",
       x = "Methylation Difference (%)",
       y = "-log10(q-value)") +
  theme_minimal()
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
DMLs <- methylKit::getMethylDiff(DMLStats_Tissue, difference = 2, qvalue = 0.05) #Identify DML based on difference threshold

length(DMLs$chr) #DML
```

    ## [1] 3

``` r
head(DMLs)
```

    ##                                       chr   start     end strand       pvalue
    ## 2122 Pocillopora_acuta_HIv2___xfSc0000006 3476631 3476631      + 3.957521e-05
    ## 2319 Pocillopora_acuta_HIv2___xfSc0000018 2653917 2653917      + 1.626684e-05
    ## 2928 Pocillopora_acuta_HIv2___xpSc0000391 1758128 1758128      + 2.866148e-05
    ##       qvalue meth.diff
    ## 2122 0.03725  41.55973
    ## 2319 0.03725  14.08451
    ## 2928 0.03725 -22.03390

## 0.4 Further look at genome wide methylation

``` r
diffMethPerChr(DMLStats_Tissue, meth.cutoff = 2, qvalue.cutoff = 0.05,cex.names=.75)
```

    ## Warning in eval(quote(list(...)), env): NAs introduced by coercion

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### 0.4.1 Annotation

``` r
require("genomationData")
require("genomation")

gff.file = "../references/Pocillopora_acuta_HIv2.gtf"
gff = gffToGRanges(gff.file)
head(gff)
```

    ## GRanges object with 6 ranges and 6 metadata columns:
    ##                     seqnames    ranges strand |   source       type     score
    ##                        <Rle> <IRanges>  <Rle> | <factor>   <factor> <numeric>
    ##   [1] Pocillopora_acuta_HI..  151-2746      + | AUGUSTUS transcript        NA
    ##   [2] Pocillopora_acuta_HI..   151-172      + | AUGUSTUS exon              NA
    ##   [3] Pocillopora_acuta_HI..   264-304      + | AUGUSTUS exon              NA
    ##   [4] Pocillopora_acuta_HI.. 1491-1602      + | AUGUSTUS exon              NA
    ##   [5] Pocillopora_acuta_HI.. 1889-1990      + | AUGUSTUS exon              NA
    ##   [6] Pocillopora_acuta_HI.. 2107-2127      + | AUGUSTUS exon              NA
    ##           phase          transcript_id                gene_id
    ##       <integer>            <character>            <character>
    ##   [1]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   [2]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   [3]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   [4]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   [5]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   [6]      <NA> Pocillopora_acuta_HI.. Pocillopora_acuta_HI..
    ##   -------
    ##   seqinfo: 425 sequences from an unspecified genome; no seqlengths

``` r
exons = gffToGRanges(gff.file, filter = "exon")

transcripts = gffToGRanges(gff.file, filter = "transcript")
```

``` r
DML_grange = as(DMLs,"GRanges")
head(DML_grange)
```

    ## GRanges object with 3 ranges and 3 metadata columns:
    ##                     seqnames    ranges strand |      pvalue    qvalue meth.diff
    ##                        <Rle> <IRanges>  <Rle> |   <numeric> <numeric> <numeric>
    ##   [1] Pocillopora_acuta_HI..   3476631      + | 3.95752e-05   0.03725   41.5597
    ##   [2] Pocillopora_acuta_HI..   2653917      + | 1.62668e-05   0.03725   14.0845
    ##   [3] Pocillopora_acuta_HI..   1758128      + | 2.86615e-05   0.03725  -22.0339
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths

``` r
transcripts = gffToGRanges(gff.file, filter = "transcript")

# Find overlaps between DMLs and transcripts
overlaps_transcripts <- findOverlaps(DML_grange, transcripts,ignore.strand = TRUE)

# Extract matching transcript information
DML_transcript_annot <- data.frame(
  DML_chr = seqnames(DML_grange)[queryHits(overlaps_transcripts)],
  DML_start = start(DML_grange)[queryHits(overlaps_transcripts)],
  DML_end = end(DML_grange)[queryHits(overlaps_transcripts)],
  DML_qvalue = (DML_grange$qvalue)[queryHits(overlaps_transcripts)],
  DML_methdiff = (DML_grange$meth.diff)[queryHits(overlaps_transcripts)],
  transcript_chr = seqnames(transcripts)[subjectHits(overlaps_transcripts)],
  transcript_start = start(transcripts)[subjectHits(overlaps_transcripts)],
  transcript_end = end(transcripts)[subjectHits(overlaps_transcripts)],
  transcript_id = transcripts$transcript_id[subjectHits(overlaps_transcripts)],
  gene_id = transcripts$gene_id[subjectHits(overlaps_transcripts)]
)

# View first few rows
head(DML_transcript_annot)
```

    ##  [1] DML_chr          DML_start        DML_end          DML_qvalue      
    ##  [5] DML_methdiff     transcript_chr   transcript_start transcript_end  
    ##  [9] transcript_id    gene_id         
    ## <0 rows> (or 0-length row.names)

## 0.5 Are any DMGs DMLs?

``` r
#load in DESeq results
DESeq <- read.csv("../output_RNA/differential_expression/DESeq_results.csv", header = TRUE) %>% dplyr::rename("query" ="X")

#make dataframes of just differentially expressed genes for each LFC direction
DE_05_Aboral <- DESeq %>% filter(padj < 0.05 & log2FoldChange > 0)
DE_05_OralEpi <- DESeq %>% filter(padj < 0.05& log2FoldChange < 0)
DE_05 <- DESeq %>% filter(padj < 0.05)

DML_transcript_annot[DML_transcript_annot$transcript_id %in% DE_05$query,]
```

    ##  [1] DML_chr          DML_start        DML_end          DML_qvalue      
    ##  [5] DML_methdiff     transcript_chr   transcript_start transcript_end  
    ##  [9] transcript_id    gene_id         
    ## <0 rows> (or 0-length row.names)

``` r
DE_05[DE_05$query %in% DML_transcript_annot$transcript_id,]
```

    ## [1] query          baseMean       log2FoldChange lfcSE          pvalue        
    ## [6] padj          
    ## <0 rows> (or 0-length row.names)

``` r
plot_data <- merge(DML_transcript_annot, DESeq, by.x = "transcript_id", by.y = "query")

ggplot(plot_data, aes(x = DML_methdiff, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), size = 3) +
  geom_text_repel(aes(label = ifelse(padj < 0.05, transcript_id, "")), max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Methylation Difference (%)", 
       y = "RNA-seq: Log2FoldChange",
       title = "Differentially methylated loci (DMLs) vs. Gene Expression",
       color = "Significant DEGs") +
  theme_minimal()
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
# Find overlaps between methylated loci and transcripts
MLs_in_genes <- regionCounts(meth_filter_destrand, regions=transcripts)
```

    ## Warning in .merge_two_Seqinfo_objects(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': Pocillopora_acuta_HIv2___Sc0000026, Pocillopora_acuta_HIv2___Sc0000028, Pocillopora_acuta_HIv2___Sc0000036, Pocillopora_acuta_HIv2___Sc0000037, Pocillopora_acuta_HIv2___Sc0000041, Pocillopora_acuta_HIv2___Sc0000046, Pocillopora_acuta_HIv2___Sc0000047, Pocillopora_acuta_HIv2___Sc0000050, Pocillopora_acuta_HIv2___Sc0000052, Pocillopora_acuta_HIv2___Sc0000053, Pocillopora_acuta_HIv2___Sc0000056, Pocillopora_acuta_HIv2___Sc0000058, Pocillopora_acuta_HIv2___Sc0000059, Pocillopora_acuta_HIv2___Sc0000060, Pocillopora_acuta_HIv2___Sc0000062, Pocillopora_acuta_HIv2___Sc0000063, Pocillopora_acuta_HIv2___Sc0000064, Pocillopora_acuta_HIv2___Sc0000066, Pocillopora_acuta_HIv2___Sc0000068, Pocillopora_acuta_HIv2___Sc0000070, Pocillopora_acuta_HIv2___Sc0000073, Pocillopora_acuta_HIv2___Sc0000075, Pocillopora_acuta_HIv2___Sc0000076, Pocillopora_acuta_HIv2___Sc0000077, Pocillopora_acuta_HIv2___Sc0000078, Pocillopora_acuta_HIv2___Sc0000080, Pocillopora_acuta_HIv2___xfSc0000031, Pocillopora_acuta_HIv2___xfSc0000032, Pocillopora_acuta_HIv2___xfSc0000034, Pocillopora_acuta_HIv2___xfSc0000041, Pocillopora_acuta_HIv2___xfSc0000042, Pocillopora_acuta_HIv2___xfSc0000044, Pocillopora_acuta_HIv2___xfSc0000045, Pocillopora_acuta_HIv2___xfSc0000047, Pocillopora_acuta_HIv2___xfSc0000050, Pocillopora_acuta_HIv2___xfSc0000051, Pocillopora_acuta_HIv2___xfSc0000055, Pocillopora_acuta_HIv2___xfSc0000056, Pocillopora_acuta_HIv2___xfSc0000057, Pocillopora_acuta_HIv2___xfSc0000058, Pocillopora_acuta_HIv2___xfSc0000062, Pocillopora_acuta_HIv2___xfSc0000063, Pocillopora_acuta_HIv2___xfSc0000064, Pocillopora_acuta_HIv2___xfSc0000065, Pocillopora_acuta_HIv2___xfSc0000066, Pocillopora_acuta_HIv2___xfSc0000070, Pocillopora_acuta_HIv2___xfSc0000071, Pocillopora_acuta_HIv2___xfSc0000072, Pocillopora_acuta_HIv2___xfSc0000073, Pocillopora_acuta_HIv2___xfSc0000074, Pocillopora_acuta_HIv2___xfSc0000075, Pocillopora_acuta_HIv2___xfSc0000077, Pocillopora_acuta_HIv2___xfSc0000078, Pocillopora_acuta_HIv2___xfSc0000079, Pocillopora_acuta_HIv2___xfSc0000080, Pocillopora_acuta_HIv2___xfSc0000081, Pocillopora_acuta_HIv2___xfSc0000084, Pocillopora_acuta_HIv2___xfSc0000085, Pocillopora_acuta_HIv2___xfSc0000086, Pocillopora_acuta_HIv2___xfSc0000087, Pocillopora_acuta_HIv2___xfSc0000089, Pocillopora_acuta_HIv2___xfSc0000090, Pocillopora_acuta_HIv2___xfSc0000091, Pocillopora_acuta_HIv2___xfSc0000092, Pocillopora_acuta_HIv2___xfSc0000094, Pocillopora_acuta_HIv2___xfSc0000095, Pocillopora_acuta_HIv2___xfSc0000096, Pocillopora_acuta_HIv2___xfSc0000097, Pocillopora_acuta_HIv2___xfSc0000098, Pocillopora_acuta_HIv2___xfSc0000099, Pocillopora_acuta_HIv2___xfSc0000100, Pocillopora_acuta_HIv2___xfSc0000101, Pocillopora_acuta_HIv2___xfSc0000102, Pocillopora_acuta_HIv2___xfSc0000103, Pocillopora_acuta_HIv2___xfSc0000104, Pocillopora_acuta_HIv2___xfSc0000105, Pocillopora_acuta_HIv2___xfSc0000107, Pocillopora_acuta_HIv2___xfSc0000108, Pocillopora_acuta_HIv2___xfSc0000109, Pocillopora_acuta_HIv2___xfSc0000110, Pocillopora_acuta_HIv2___xfSc0000113, Pocillopora_acuta_HIv2___xfSc0000114, Pocillopora_acuta_HIv2___xfSc0000115, Pocillopora_acuta_HIv2___xfSc0000116, Pocillopora_acuta_HIv2___xfSc0000117, Pocillopora_acuta_HIv2___xfSc0000119, Pocillopora_acuta_HIv2___xfSc0000121, Pocillopora_acuta_HIv2___xfSc0000122, Pocillopora_acuta_HIv2___xfSc0000125, Pocillopora_acuta_HIv2___xfSc0000126, Pocillopora_acuta_HIv2___xfSc0000127, Pocillopora_acuta_HIv2___xfSc0000128, Pocillopora_acuta_HIv2___xfSc0000129, Pocillopora_acuta_HIv2___xfSc0000130, Pocillopora_acuta_HIv2___xfSc0000131, Pocillopora_acuta_HIv2___xfSc0000132, Pocillopora_acuta_HIv2___xfSc0000133, Pocillopora_acuta_HIv2___xfSc0000135, Pocillopora_acuta_HIv2___xfSc0000136, Pocillopora_acuta_HIv2___xfSc0000137, Pocillopora_acuta_HIv2___xfSc0000138, Pocillopora_acuta_HIv2___xfSc0000139, Pocillopora_acuta_HIv2___xfSc0000140, Pocillopora_acuta_HIv2___xfSc0000141, Pocillopora_acuta_HIv2___xfSc0000142, Pocillopora_acuta_HIv2___xfSc0000143, Pocillopora_acuta_HIv2___xfSc0000144, Pocillopora_acuta_HIv2___xfSc0000145, Pocillopora_acuta_HIv2___xfSc0000147, Pocillopora_acuta_HIv2___xfSc0000148, Pocillopora_acuta_HIv2___xfSc0000150, Pocillopora_acuta_HIv2___xfSc0000151, Pocillopora_acuta_HIv2___xfSc0000152, Pocillopora_acuta_HIv2___xfSc0000154, Pocillopora_acuta_HIv2___xfSc0000155, Pocillopora_acuta_HIv2___xfSc0000157, Pocillopora_acuta_HIv2___xfSc0000158, Pocillopora_acuta_HIv2___xfSc0000159, Pocillopora_acuta_HIv2___xfSc0000160, Pocillopora_acuta_HIv2___xfSc0000161, Pocillopora_acuta_HIv2___xfSc0000162, Pocillopora_acuta_HIv2___xfSc0000163, Pocillopora_acuta_HIv2___xfSc0000164, Pocillopora_acuta_HIv2___xfSc0000165, Pocillopora_acuta_HIv2___xfSc0000166, Pocillopora_acuta_HIv2___xfSc0000167, Pocillopora_acuta_HIv2___xfSc0000168, Pocillopora_acuta_HIv2___xfSc0000170, Pocillopora_acuta_HIv2___xfSc0000171, Pocillopora_acuta_HIv2___xfSc0000173, Pocillopora_acuta_HIv2___xfSc0000174, Pocillopora_acuta_HIv2___xfSc0000177, Pocillopora_acuta_HIv2___xfSc0000178, Pocillopora_acuta_HIv2___xfSc0000179, Pocillopora_acuta_HIv2___xfSc0000180, Pocillopora_acuta_HIv2___xfSc0000181, Pocillopora_acuta_HIv2___xfSc0000182, Pocillopora_acuta_HIv2___xfSc0000183, Pocillopora_acuta_HIv2___xfSc0000184, Pocillopora_acuta_HIv2___xfSc0000185, Pocillopora_acuta_HIv2___xfSc0000186, Pocillopora_acuta_HIv2___xfSc0000187, Pocillopora_acuta_HIv2___xfSc0000188, Pocillopora_acuta_HIv2___xfSc0000189, Pocillopora_acuta_HIv2___xfSc0000190, Pocillopora_acuta_HIv2___xfSc0000191, Pocillopora_acuta_HIv2___xfSc0000192, Pocillopora_acuta_HIv2___xfSc0000193, Pocillopora_acuta_HIv2___xfSc0000195, Pocillopora_acuta_HIv2___xfSc0000196, Pocillopora_acuta_HIv2___xfSc0000197, Pocillopora_acuta_HIv2___xfSc0000200, Pocillopora_acuta_HIv2___xfSc0000202, Pocillopora_acuta_HIv2___xfSc0000203, Pocillopora_acuta_HIv2___xfSc0000204, Pocillopora_acuta_HIv2___xfSc0000205, Pocillopora_acuta_HIv2___xfSc0000206, Pocillopora_acuta_HIv2___xfSc0000207, Pocillopora_acuta_HIv2___xfSc0000208, Pocillopora_acuta_HIv2___xfSc0000209, Pocillopora_acuta_HIv2___xfSc0000210, Pocillopora_acuta_HIv2___xfSc0000212, Pocillopora_acuta_HIv2___xfSc0000213, Pocillopora_acuta_HIv2___xfSc0000214, Pocillopora_acuta_HIv2___xfSc0000215, Pocillopora_acuta_HIv2___xfSc0000216, Pocillopora_acuta_HIv2___xfSc0000217, Pocillopora_acuta_HIv2___xfSc0000218, Pocillopora_acuta_HIv2___xfSc0000220, Pocillopora_acuta_HIv2___xfSc0000221, Pocillopora_acuta_HIv2___xfSc0000222, Pocillopora_acuta_HIv2___xfSc0000223, Pocillopora_acuta_HIv2___xfSc0000224, Pocillopora_acuta_HIv2___xfSc0000226, Pocillopora_acuta_HIv2___xfSc0000227, Pocillopora_acuta_HIv2___xfSc0000228, Pocillopora_acuta_HIv2___xfSc0000230, Pocillopora_acuta_HIv2___xfSc0000231, Pocillopora_acuta_HIv2___xfSc0000232, Pocillopora_acuta_HIv2___xfSc0000233, Pocillopora_acuta_HIv2___xfSc0000234, Pocillopora_acuta_HIv2___xfSc0000235, Pocillopora_acuta_HIv2___xfSc0000236, Pocillopora_acuta_HIv2___xfSc0000237, Pocillopora_acuta_HIv2___xfSc0000238, Pocillopora_acuta_HIv2___xfSc0000239, Pocillopora_acuta_HIv2___xfSc0000240, Pocillopora_acuta_HIv2___xfSc0000241, Pocillopora_acuta_HIv2___xfSc0000242, Pocillopora_acuta_HIv2___xfSc0000243, Pocillopora_acuta_HIv2___xfSc0000244, Pocillopora_acuta_HIv2___xfSc0000245, Pocillopora_acuta_HIv2___xfSc0000247, Pocillopora_acuta_HIv2___xfSc0000248, Pocillopora_acuta_HIv2___xfSc0000249, Pocillopora_acuta_HIv2___xfSc0000250, Pocillopora_acuta_HIv2___xfSc0000251, Pocillopora_acuta_HIv2___xfSc0000252, Pocillopora_acuta_HIv2___xfSc0000254, Pocillopora_acuta_HIv2___xfSc0000255, Pocillopora_acuta_HIv2___xfSc0000256, Pocillopora_acuta_HIv2___xfSc0000257, Pocillopora_acuta_HIv2___xfSc0000258, Pocillopora_acuta_HIv2___xfSc0000259, Pocillopora_acuta_HIv2___xfSc0000260, Pocillopora_acuta_HIv2___xfSc0000261, Pocillopora_acuta_HIv2___xfSc0000262, Pocillopora_acuta_HIv2___xfSc0000264, Pocillopora_acuta_HIv2___xfSc0000265, Pocillopora_acuta_HIv2___xfSc0000266, Pocillopora_acuta_HIv2___xfSc0000267, Pocillopora_acuta_HIv2___xfSc0000268, Pocillopora_acuta_HIv2___xfSc0000269, Pocillopora_acuta_HIv2___xfSc0000270, Pocillopora_acuta_HIv2___xfSc

``` r
#MLs_in_genes <- getData(MLs_in_genes)
MLs_in_genes$gene_id <- transcripts$gene_id[match(paste(MLs_in_genes$chr, MLs_in_genes$start, MLs_in_genes$end), paste(seqnames(transcripts), start(transcripts), end(transcripts)))]

percent_meth <- percMethylation(MLs_in_genes)
percent_meth <- as.data.frame(percent_meth)
percent_meth$gene_id <- MLs_in_genes$gene_id 

percent_meth <- percent_meth %>% rowwise() %>%
  mutate(percent_meth_ALL = mean(c_across(starts_with("LCM")), na.rm = TRUE)) %>%
  ungroup()

percent_meth <- percent_meth %>% rowwise() %>%
  mutate(Oral = mean(c_across(meta$Sample[meta$Tissue=="OralEpi"]), na.rm = TRUE)) %>%
  mutate(Aboral = mean(c_across(meta$Sample[meta$Tissue=="Aboral"]), na.rm = TRUE)) %>%
  ungroup()
  
percent_meth_long <- percent_meth %>% pivot_longer(
                                      cols = c(Oral,Aboral),
                                      names_to = "Tissue",
                                      values_to = "tissue_percent_meth"
)
```

``` r
plot_data <- merge(percent_meth, DESeq, by.x = "gene_id", by.y = "query")

ggplot(plot_data, aes(y = baseMean, x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "DeSeq2 BaseMean", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ggplot(plot_data, aes(y = log2(baseMean), x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Log2(DeSeq2 BaseMean)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

``` r
ggplot(plot_data, aes(y = abs(log2FoldChange), x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +  stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Abs(DeSeq2 Log2FoldChange)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->

``` r
# Create the plot
ggplot(plot_data, aes(x = percent_meth_ALL, y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05)) +
  #geom_text_repel(aes(label = ifelse(padj < 0.05, gene_id, "")), max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Average CpG % methylation of gene", 
       y = "RNA-seq: Log2FoldChange",
       title = "Differentially methylated loci (DMLs) vs. Gene Expression",
       color = "Significant DEGs") +
  theme_minimal()
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-21-4.png)<!-- -->

``` r
plot_data <- merge(percent_meth_long, DESeq, by.x = "gene_id", by.y = "query")

ggplot(plot_data, aes(y = abs(log2FoldChange), x = tissue_percent_meth, color=Tissue)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2"))+
  labs(x = "Average CpG % methylation of gene", y = "Abs(DeSeq2 Log2FoldChange)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### 0.5.1 Read counts

RNA metadata:

``` r
meta_RNA <- read.csv("../data_RNA/LCM_RNA_metadata.csv") %>%
            dplyr::arrange(Sample) %>%
            mutate(across(c(Tissue, Fragment, Section_Date, LCM_Date), factor)) # Set variables as factors 

meta_RNA$Tissue <- factor(meta_RNA$Tissue, levels = c("OralEpi","Aboral")) #we want OralEpi to be the baseline
```

``` r
filtered_counts_RNA <- read.csv("../output_RNA/differential_expression/filtered_counts.csv") 
rownames(filtered_counts_RNA) <- filtered_counts_RNA$X
filtered_counts_RNA <- filtered_counts_RNA %>% dplyr::select(-X)
filtered_counts_RNA$sum <- rowSums(filtered_counts_RNA)
filtered_counts_RNA$gene_id <- rownames(filtered_counts_RNA)

filtered_counts_RNA <- filtered_counts_RNA %>% rowwise() %>%
  mutate(Oral = mean(c_across(meta_RNA$Sample[meta_RNA$Tissue=="OralEpi"]), na.rm = TRUE)) %>%
  mutate(Aboral = mean(c_across(meta_RNA$Sample[meta_RNA$Tissue=="Aboral"]), na.rm = TRUE)) %>%
  ungroup()
  
filtered_counts_RNA_long <- filtered_counts_RNA %>% pivot_longer(
                                      cols = c(Oral,Aboral),
                                      names_to = "Tissue",
                                      values_to = "tissue_mean_counts"
)
```

``` r
plot_data <- merge(percent_meth_long, filtered_counts_RNA_long, by = c("gene_id", "Tissue"))

ggplot(plot_data, aes(y = log2(tissue_mean_counts), x = tissue_percent_meth, color=Tissue)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2"))+
  labs(x = "Average CpG % methylation of gene", y = "Mean counts, all samples", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 3 rows containing non-finite outside the scale range
    ## (`stat_poly_eq()`).

    ## Warning in ci_f_ncp(stat, df1 = df1, df2 = df2, probs = probs): Upper limit
    ## outside search range. Set to the maximum of the parameter range.

    ## Warning in compute_group(...): CI computation error: Error in check_output(cint, probs = probs, parameter_range = c(0, 1)): out[1] <= out[2] is not TRUE

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->
