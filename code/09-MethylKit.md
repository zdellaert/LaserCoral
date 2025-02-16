09-MethylKit
================
Zoe Dellaert
2025-02-09

- [0.1 MethylKit](#01-methylkit)
- [0.2 Managing Packages Using Renv](#02-managing-packages-using-renv)
- [0.3 Load packages](#03-load-packages)
  - [0.3.1 other possible filtering](#031-other-possible-filtering)
  - [0.3.2 Identify DML](#032-identify-dml)
- [0.4 Further look at genome wide
  methylation](#04-further-look-at-genome-wide-methylation)
  - [0.4.1 Annotation](#041-annotation)
- [0.5 Are any DMGs DMLs?](#05-are-any-dmgs-dmls)

## 0.1 MethylKit

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
    ##  [4] gplots_3.2.0          vegan_2.6-10          lattice_0.22-6       
    ##  [7] permute_0.9-7         lubridate_1.9.4       forcats_1.0.0        
    ## [10] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.4          
    ## [13] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
    ## [16] ggplot2_3.5.1         tidyverse_2.0.0       methylKit_1.32.0     
    ## [19] GenomicRanges_1.58.0  GenomeInfoDb_1.42.3   IRanges_2.40.1       
    ## [22] S4Vectors_0.44.0      BiocGenerics_0.52.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-9                rlang_1.1.5                
    ##  [3] magrittr_2.0.3              gridBase_0.4-7             
    ##  [5] matrixStats_1.5.0           compiler_4.4.0             
    ##  [7] mgcv_1.9-1                  vctrs_0.6.5                
    ##  [9] reshape2_1.4.4              pkgconfig_2.0.3            
    ## [11] crayon_1.5.3                fastmap_1.2.0              
    ## [13] XVector_0.46.0              caTools_1.18.3             
    ## [15] Rsamtools_2.22.0            rmarkdown_2.29             
    ## [17] tzdb_0.4.0                  UCSC.utils_1.2.0           
    ## [19] xfun_0.50                   zlibbioc_1.52.0            
    ## [21] jsonlite_1.8.9              DelayedArray_0.32.0        
    ## [23] BiocParallel_1.40.0         parallel_4.4.0             
    ## [25] cluster_2.1.8               R6_2.6.0                   
    ## [27] stringi_1.8.4               limma_3.62.2               
    ## [29] rtracklayer_1.66.0          numDeriv_2016.8-1.1        
    ## [31] Rcpp_1.0.14                 SummarizedExperiment_1.36.0
    ## [33] knitr_1.49                  R.utils_2.12.3             
    ## [35] Matrix_1.7-2                splines_4.4.0              
    ## [37] timechange_0.3.0            tidyselect_1.2.1           
    ## [39] seqPattern_1.38.0           qvalue_2.38.0              
    ## [41] rstudioapi_0.17.1           abind_1.4-8                
    ## [43] yaml_2.3.10                 codetools_0.2-20           
    ## [45] curl_6.2.0                  plyr_1.8.9                 
    ## [47] Biobase_2.66.0              withr_3.0.2                
    ## [49] coda_0.19-4.1               evaluate_1.0.3             
    ## [51] mclust_6.1.1                Biostrings_2.74.1          
    ## [53] pillar_1.10.1               BiocManager_1.30.25        
    ## [55] MatrixGenerics_1.18.1       KernSmooth_2.23-26         
    ## [57] renv_1.1.1                  generics_0.1.3             
    ## [59] RCurl_1.98-1.16             emdbook_1.3.13             
    ## [61] hms_1.1.3                   munsell_0.5.1              
    ## [63] scales_1.3.0                gtools_3.9.5               
    ## [65] glue_1.8.0                  tools_4.4.0                
    ## [67] BiocIO_1.16.0               data.table_1.16.4          
    ## [69] BSgenome_1.74.0             GenomicAlignments_1.42.0   
    ## [71] mvtnorm_1.3-3               XML_3.99-0.18              
    ## [73] plotrix_3.8-4               impute_1.80.0              
    ## [75] bbmle_1.0.25.1              bdsmatrix_1.3-7            
    ## [77] colorspace_2.1-1            nlme_3.1-167               
    ## [79] GenomeInfoDbData_1.2.13     restfulr_0.0.15            
    ## [81] cli_3.6.4                   fastseg_1.52.0             
    ## [83] S4Arrays_1.6.0              gtable_0.3.6               
    ## [85] R.methodsS3_1.8.2           digest_0.6.37              
    ## [87] SparseArray_1.6.1           rjson_0.2.23               
    ## [89] htmltools_0.5.8.1           R.oo_1.27.0                
    ## [91] lifecycle_1.0.4             httr_1.4.7                 
    ## [93] statmod_1.5.0               MASS_7.3-64

``` r
meta <- read.csv("../data_WGBS/LCM_WGBS_metadata.csv", sep = ",", header = TRUE) %>%
  mutate(Section_Date = as.character(Section_Date), LCM_Date = as.character(LCM_Date),DNA_Extraction_Date = as.character(DNA_Extraction_Date))

meta <- meta %>% arrange(Sample)
```

``` r
file_list <- list.files("../output_WGBS/methylseq_V3_bwa_test/methyldackel", pattern = "*methylKit", full.names = TRUE, include.dirs = FALSE)

sample <- gsub(".markdup.sorted_CpG.methylKit", "", basename(file_list) )

sample == meta$Sample #the files and metadata are in the same order
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
tissue <- meta$Tissue
tissue_binary <- gsub("Aboral", "1", tissue)
tissue_binary <- gsub("OralEpi", "0", tissue_binary)
tissue_binary <- as.numeric(tissue_binary)
fragment <- meta$Fragment

# methylObj=methRead(as.list(file_list),
#            sample.id = as.list(sample),
#            assembly = "Pacuta",
#            treatment = tissue_binary,
#            context = "CpG",
#            mincov = 10
#            )
#   
# save(methylObj, file = "../output_WGBS/MethylKit.RData") 
```

``` r
load("../output_WGBS/MethylKit.RData")

getMethylationStats(methylObj[[2]],plot=FALSE,both.strands=FALSE)
```

    ## methylation statistics per base
    ## summary:
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.00    0.00    0.00   12.70   18.18  100.00 
    ## percentiles:
    ##        0%       10%       20%       30%       40%       50%       60%       70% 
    ##   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000   9.52381 
    ##       80%       90%       95%       99%     99.5%     99.9%      100% 
    ##  27.27273  45.45455  64.28571 100.00000 100.00000 100.00000 100.00000

``` r
getMethylationStats(methylObj[[2]],plot=TRUE,both.strands=FALSE)
```

![](09-MethylKit_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
getCoverageStats(methylObj[[2]],plot=TRUE,both.strands=FALSE)
```

![](09-MethylKit_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
filtered_methylObj=filterByCoverage(methylObj,lo.count=10,lo.perc=NULL,
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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
PCASamples(meth_filter_destrand)
```

![](09-MethylKit_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
getCorrelation(meth_filter_destrand,plot=TRUE)
```

    ##            LCM_1    LCM_11    LCM_12    LCM_17    LCM_18    LCM_24    LCM_25
    ## LCM_1  1.0000000 0.8071167 0.7926500 0.7987458 0.7720128 0.6870701 0.7760903
    ## LCM_11 0.8071167 1.0000000 0.8090847 0.8112960 0.7812533 0.6885644 0.7713116
    ## LCM_12 0.7926500 0.8090847 1.0000000 0.8166968 0.7694800 0.6984713 0.7666580
    ## LCM_17 0.7987458 0.8112960 0.8166968 1.0000000 0.7715643 0.6852749 0.7733168
    ## LCM_18 0.7720128 0.7812533 0.7694800 0.7715643 1.0000000 0.7011845 0.7711137
    ## LCM_24 0.6870701 0.6885644 0.6984713 0.6852749 0.7011845 1.0000000 0.7222715
    ## LCM_25 0.7760903 0.7713116 0.7666580 0.7733168 0.7711137 0.7222715 1.0000000
    ## LCM_3  0.7734702 0.7998513 0.7887419 0.7740742 0.7522157 0.6579610 0.7343776
    ## LCM_32 0.7591867 0.7810416 0.7644149 0.7624203 0.7641989 0.6607244 0.7409674
    ## LCM_33 0.7685862 0.8026637 0.7790518 0.7669838 0.7445778 0.6355911 0.7290590
    ##            LCM_3    LCM_32    LCM_33
    ## LCM_1  0.7734702 0.7591867 0.7685862
    ## LCM_11 0.7998513 0.7810416 0.8026637
    ## LCM_12 0.7887419 0.7644149 0.7790518
    ## LCM_17 0.7740742 0.7624203 0.7669838
    ## LCM_18 0.7522157 0.7641989 0.7445778
    ## LCM_24 0.6579610 0.6607244 0.6355911
    ## LCM_25 0.7343776 0.7409674 0.7290590
    ## LCM_3  1.0000000 0.7600768 0.7669578
    ## LCM_32 0.7600768 1.0000000 0.7727821
    ## LCM_33 0.7669578 0.7727821 1.0000000

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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-7-1.png)<!-- --> \###
batch effects

``` r
as=assocComp(mBase=meth_filter_destrand,select(meta,c("PCR_ReAmp_Cycles", "Fragment")))
as
```

    ## $pcs
    ##               PC1         PC2         PC3         PC4          PC5         PC6
    ## LCM_1  -0.3216779  0.05915356  0.30196916 -0.14420962  0.106034156  0.02317996
    ## LCM_11 -0.3265585  0.15514843  0.09086618  0.10898681  0.060885171 -0.09416565
    ## LCM_12 -0.3237446  0.07166955  0.28814038  0.19168049  0.002164678 -0.12638941
    ## LCM_17 -0.3227715  0.08032452  0.38387620 -0.05271814  0.062569105 -0.12800268
    ## LCM_18 -0.3170704 -0.09337650 -0.24292276 -0.48602692 -0.476847742 -0.58246573
    ## LCM_24 -0.2877121 -0.83279369 -0.15221492  0.42336874  0.065202942 -0.06486596
    ## LCM_25 -0.3152164 -0.25590050  0.16930676 -0.54083301  0.234757957  0.51179484
    ## LCM_3  -0.3164516  0.22392814  0.07429410  0.39346280 -0.627734867  0.38598279
    ## LCM_32 -0.3146195  0.17257995 -0.69944046 -0.09631807  0.031696993  0.35712452
    ## LCM_33 -0.3147866  0.33639865 -0.25414553  0.23611671  0.547108693 -0.27623910
    ##                 PC7         PC8         PC9        PC10
    ## LCM_1   0.139546895  0.82219692  0.26173767 -0.07262323
    ## LCM_11  0.109075131  0.06607662 -0.56734877  0.70599677
    ## LCM_12 -0.367549042 -0.30982410  0.63482906  0.34672241
    ## LCM_17 -0.521088871 -0.10202600 -0.43947948 -0.49380329
    ## LCM_18  0.139134965 -0.07779117  0.05051128 -0.03003498
    ## LCM_24  0.001981578  0.08167280 -0.05123066 -0.05922844
    ## LCM_25  0.263831482 -0.35517125  0.01959832  0.02837934
    ## LCM_3   0.305140224 -0.09386892 -0.01449674 -0.21027636
    ## LCM_32 -0.463220542  0.16450265  0.02599882  0.04340127
    ## LCM_33  0.406632547 -0.19617546  0.08374050 -0.28446617
    ## 
    ## $vars
    ##  [1] 78.100976  4.266914  2.673250  2.621410  2.334903  2.187428  2.181249
    ##  [8]  2.072173  1.829281  1.732416
    ## 
    ## $association
    ##                         PC1       PC2       PC3       PC4       PC5       PC6
    ## PCR_ReAmp_Cycles 0.78459351 0.3339619 0.1283833 0.4572514 0.7734739 0.2419045
    ## Fragment         0.09291936 0.1553095 0.3205067 0.8780986 0.5487392 0.2249586
    ##                        PC7       PC8       PC9      PC10
    ## PCR_ReAmp_Cycles 0.4221351 0.7070921 0.6163093 0.1533615
    ## Fragment         0.5855722 0.8410289 0.8410289 0.2536518

### 0.3.1 other possible filtering

``` r
# get percent methylation matrix
pm=percMethylation(meth_filter_destrand)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)
```

![](09-MethylKit_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# keep only CpG with standard deviations larger than 2%
meth <- meth_filter_destrand[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth_filter_destrand)
```

    ## [1] 5669

``` r
nrow(meth)
```

    ## [1] 5481

### 0.3.2 Identify DML

``` r
DMLStats_Tissue <- methylKit::calculateDiffMeth(meth_filter_destrand, overdispersion = "MN", test = "Chisq", mc.cores = 8) #Calculate differential methylation statistics and include covariate information.
```

    ## two groups detected:
    ##  will calculate methylation difference as the difference of
    ## treatment (group: 1) - control (group: 0)

``` r
head(DMLStats_Tissue) #Look at differential methylation output
```

    ##                                  chr  start    end strand     pvalue    qvalue
    ## 1 Pocillopora_acuta_HIv2___Sc0000000  48656  48656      + 0.34173181 0.9030236
    ## 2 Pocillopora_acuta_HIv2___Sc0000000 508417 508417      + 0.06541963 0.6828262
    ## 3 Pocillopora_acuta_HIv2___Sc0000000 508643 508643      + 0.17842512 0.8308404
    ## 4 Pocillopora_acuta_HIv2___Sc0000000 509163 509163      + 0.30851610 0.8891573
    ## 5 Pocillopora_acuta_HIv2___Sc0000000 844614 844614      + 0.75207902 0.9655994
    ## 6 Pocillopora_acuta_HIv2___Sc0000000 844624 844624      + 0.62496289 0.9494677
    ##   meth.diff
    ## 1  7.382703
    ## 2 19.510490
    ## 3  8.511694
    ## 4  6.805556
    ## 5 -2.404868
    ## 6 -3.431472

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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

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

![](09-MethylKit_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
DMLs <- methylKit::getMethylDiff(DMLStats_Tissue, difference = 2, qvalue = 0.05) #Identify DML based on difference threshold

length(DMLs$chr) #DML
```

    ## [1] 9

``` r
head(DMLs)
```

    ##                                     chr   start     end strand       pvalue
    ## 163  Pocillopora_acuta_HIv2___Sc0000001 3408801 3408801      + 2.890498e-05
    ## 166  Pocillopora_acuta_HIv2___Sc0000001 3408829 3408829      + 3.437585e-06
    ## 170  Pocillopora_acuta_HIv2___Sc0000001 3408857 3408857      + 1.639072e-06
    ## 670  Pocillopora_acuta_HIv2___Sc0000003 9609977 9609977      + 5.074907e-05
    ## 2203 Pocillopora_acuta_HIv2___Sc0000020 3271786 3271786      + 1.874037e-05
    ## 2448 Pocillopora_acuta_HIv2___Sc0000024  987115  987115      + 1.139080e-06
    ##           qvalue  meth.diff
    ## 163  0.022750677 -17.780636
    ## 166  0.004734927 -20.282751
    ## 170  0.003010208 -22.591117
    ## 670  0.031067429   4.918137
    ## 2203 0.017208641  18.014797
    ## 2448 0.003010208 -25.779357

## 0.4 Further look at genome wide methylation

``` r
diffMethPerChr(DMLStats_Tissue, meth.cutoff = 2, qvalue.cutoff = 0.05,cex.names=.75)
```

    ## Warning in eval(quote(list(...)), env): NAs introduced by coercion

![](09-MethylKit_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

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

    ## GRanges object with 6 ranges and 3 metadata columns:
    ##                     seqnames    ranges strand |      pvalue     qvalue
    ##                        <Rle> <IRanges>  <Rle> |   <numeric>  <numeric>
    ##   [1] Pocillopora_acuta_HI..   3408801      + | 2.89050e-05 0.02275068
    ##   [2] Pocillopora_acuta_HI..   3408829      + | 3.43759e-06 0.00473493
    ##   [3] Pocillopora_acuta_HI..   3408857      + | 1.63907e-06 0.00301021
    ##   [4] Pocillopora_acuta_HI..   9609977      + | 5.07491e-05 0.03106743
    ##   [5] Pocillopora_acuta_HI..   3271786      + | 1.87404e-05 0.01720864
    ##   [6] Pocillopora_acuta_HI..    987115      + | 1.13908e-06 0.00301021
    ##       meth.diff
    ##       <numeric>
    ##   [1] -17.78064
    ##   [2] -20.28275
    ##   [3] -22.59112
    ##   [4]   4.91814
    ##   [5]  18.01480
    ##   [6] -25.77936
    ##   -------
    ##   seqinfo: 6 sequences from an unspecified genome; no seqlengths

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

    ##                              DML_chr DML_start DML_end  DML_qvalue DML_methdiff
    ## 1 Pocillopora_acuta_HIv2___Sc0000001   3408801 3408801 0.022750677   -17.780636
    ## 2 Pocillopora_acuta_HIv2___Sc0000001   3408829 3408829 0.004734927   -20.282751
    ## 3 Pocillopora_acuta_HIv2___Sc0000001   3408857 3408857 0.003010208   -22.591117
    ## 4 Pocillopora_acuta_HIv2___Sc0000003   9609977 9609977 0.031067429     4.918137
    ## 5 Pocillopora_acuta_HIv2___Sc0000020   3271786 3271786 0.017208641    18.014797
    ## 6 Pocillopora_acuta_HIv2___Sc0000024    987115  987115 0.003010208   -25.779357
    ##                       transcript_chr transcript_start transcript_end
    ## 1 Pocillopora_acuta_HIv2___Sc0000001          3407978        3412324
    ## 2 Pocillopora_acuta_HIv2___Sc0000001          3407978        3412324
    ## 3 Pocillopora_acuta_HIv2___Sc0000001          3407978        3412324
    ## 4 Pocillopora_acuta_HIv2___Sc0000003          9605569        9649727
    ## 5 Pocillopora_acuta_HIv2___Sc0000020          3263869        3273843
    ## 6 Pocillopora_acuta_HIv2___Sc0000024           985385         998579
    ##                               transcript_id
    ## 1 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 2 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 3 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 4  Pocillopora_acuta_HIv2___RNAseq.g5543.t2
    ## 5 Pocillopora_acuta_HIv2___RNAseq.g14913.t1
    ## 6  Pocillopora_acuta_HIv2___RNAseq.g3700.t1
    ##                                     gene_id
    ## 1 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 2 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 3 Pocillopora_acuta_HIv2___RNAseq.g10458.t1
    ## 4  Pocillopora_acuta_HIv2___RNAseq.g5543.t2
    ## 5 Pocillopora_acuta_HIv2___RNAseq.g14913.t1
    ## 6  Pocillopora_acuta_HIv2___RNAseq.g3700.t1

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
