09-MethylKit, filtered for \>90% conversion efficiency
================
Zoe Dellaert
2025-02-09

- [0.1 MethylKit - Reads filtered for \>90% conversion
  efficiency](#01-methylkit---reads-filtered-for-90-conversion-efficiency)
- [0.2 Managing Packages Using Renv](#02-managing-packages-using-renv)
- [0.3 Load packages](#03-load-packages)
  - [0.3.1 Note: I changed the code below so that a CpG does not have to
    have 5X coverage in all samples to be analyzed, and only needs 5X
    coverage in 2 samples per group to be retained. This way we don’t
    remove CpGs that happen to have lower coverage in a few samples or
    one
    tissue.](#031-note-i-changed-the-code-below-so-that-a-cpg-does-not-have-to-have-5x-coverage-in-all-samples-to-be-analyzed-and-only-needs-5x-coverage-in-2-samples-per-group-to-be-retained-this-way-we-dont-remove-cpgs-that-happen-to-have-lower-coverage-in-a-few-samples-or-one-tissue)
  - [0.3.2 batch effects](#032-batch-effects)
  - [0.3.3 other possible filtering](#033-other-possible-filtering)
  - [0.3.4 Identify DML](#034-identify-dml)
- [0.4 Further look at genome wide
  methylation](#04-further-look-at-genome-wide-methylation)
  - [0.4.1 Annotation](#041-annotation)
- [0.5 Are any DMGs DMLs?](#05-are-any-dmgs-dmls)
  - [0.5.1 BaseMean: Raw gene expression levels (skewed by highly
    expressed
    genes)](#051-basemean-raw-gene-expression-levels-skewed-by-highly-expressed-genes)
  - [0.5.2 Log2 of BaseMean: Transformed gene expression
    levels](#052-log2-of-basemean-transformed-gene-expression-levels)
  - [0.5.3 Absolute value of Log2FoldChange: Relationship between
    overall methylation of a gene and whether or not it is
    differentially
    expressed](#053-absolute-value-of-log2foldchange-relationship-between-overall-methylation-of-a-gene-and-whether-or-not-it-is-differentially-expressed)
  - [0.5.4 Read counts/vsd transformed
    counts](#054-read-countsvsd-transformed-counts)

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
require("parallel")
```

    ## Loading required package: parallel

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
    ##  [1] grid      parallel  stats4    stats     graphics  grDevices datasets 
    ##  [8] utils     methods   base     
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
    ## [91] XML_3.99-0.18               MatrixModels_0.5-3         
    ## [93] mclust_6.1.1                bitops_1.0-9               
    ## [95] mvtnorm_1.3-3               scales_1.3.0               
    ## [97] crayon_1.5.3                rlang_1.1.5

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
```

``` r
# methylObj=methRead(as.list(file_list),
#            sample.id = as.list(sample),
#            assembly = "Pacuta",
#            treatment = tissue_binary,
#            context = "CpG",
#            mincov = 5
#            )

# save(methylObj, file = "../output_WGBS/MethylKit.RData")

methyl_list <- mclapply(seq_along(file_list), function(i) {
  methRead(
    file_list[[i]],
    sample.id = sample[[i]],
    assembly = "Pacuta",
    treatment = tissue_binary[[i]],
    context = "CpG",
    mincov = 5
  )
}, mc.cores = 4) 

methylObj <- new("methylRawList", methyl_list, treatment = tissue_binary)

save(methylObj, file = "../output_WGBS/MethylKit_20250513.RData")
```

``` r
#load("../output_WGBS/MethylKit_20250513.RData")

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
getCoverageStats(methylObj[[2]],plot=TRUE,both.strands=FALSE)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
getCoverageStats(methylObj[[2]],plot=TRUE,both.strands=TRUE)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

### 0.3.1 Note: I changed the code below so that a CpG does not have to have 5X coverage in all samples to be analyzed, and only needs 5X coverage in 2 samples per group to be retained. This way we don’t remove CpGs that happen to have lower coverage in a few samples or one tissue.

``` r
#filtered_methylObj=filterByCoverage(methylObj,lo.count=5,lo.perc=NULL,
#                                      hi.count=NULL,hi.perc=99.9)

#filtered_methylObj_norm <- filtered_methylObj %>% methylKit::normalizeCoverage(.)
methylObj_norm <- methylObj %>% methylKit::normalizeCoverage(.)
```

``` r
#meth_filter=methylKit::unite(filtered_methylObj_norm)
meth_filter=methylKit::unite(methylObj_norm, min.per.group = c(3L,3L))
```

    ## uniting...

    ## Warning in rowSums(ldat) >= min.per.group: longer object length is not a
    ## multiple of shorter object length

``` r
#meth_filter_destrand=methylKit::unite(filtered_methylObj_norm,destrand = TRUE)
meth_filter_destrand=methylKit::unite(methylObj_norm, min.per.group = c(3L,3L), destrand = TRUE)
```

    ## destranding...
    ## uniting...

    ## Warning in rowSums(ldat) >= min.per.group: longer object length is not a
    ## multiple of shorter object length
    ## Warning in rowSums(ldat) >= min.per.group: longer object length is not a
    ## multiple of shorter object length

``` r
clusterSamples(meth_filter, dist="correlation", method="ward", plot=TRUE)
```

    ## The "ward" method has been renamed to "ward.D"; note new "ward.D2"

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
PCASamples(meth_filter_destrand)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
getCorrelation(meth_filter_destrand,plot=TRUE)
```

    ##            LCM_1    LCM_11    LCM_12    LCM_17    LCM_18    LCM_24    LCM_25
    ## LCM_1  1.0000000 0.7222808 0.7057141 0.7223606 0.7507291 0.6701742 0.7235446
    ## LCM_11 0.7222808 1.0000000 0.6855233 0.6821552 0.6858891 0.6332845 0.6866522
    ## LCM_12 0.7057141 0.6855233 1.0000000 0.7235759 0.7371029 0.7103337 0.7362472
    ## LCM_17 0.7223606 0.6821552 0.7235759 1.0000000 0.7327728 0.6518541 0.7299340
    ## LCM_18 0.7507291 0.6858891 0.7371029 0.7327728 1.0000000 0.6837319 0.7504358
    ## LCM_24 0.6701742 0.6332845 0.7103337 0.6518541 0.6837319 1.0000000 0.7266801
    ## LCM_25 0.7235446 0.6866522 0.7362472 0.7299340 0.7504358 0.7266801 1.0000000
    ## LCM_3  0.6447084 0.6611449 0.6455462 0.6312393 0.6349158 0.6784281 0.6782781
    ## LCM_32 0.6620827 0.6091002 0.6556295 0.6316165 0.6513680 0.6415957 0.6665304
    ## LCM_33 0.6619763 0.5983036 0.6270250 0.6429158 0.6617203 0.5945263 0.6355318
    ##            LCM_3    LCM_32    LCM_33
    ## LCM_1  0.6447084 0.6620827 0.6619763
    ## LCM_11 0.6611449 0.6091002 0.5983036
    ## LCM_12 0.6455462 0.6556295 0.6270250
    ## LCM_17 0.6312393 0.6316165 0.6429158
    ## LCM_18 0.6349158 0.6513680 0.6617203
    ## LCM_24 0.6784281 0.6415957 0.5945263
    ## LCM_25 0.6782781 0.6665304 0.6355318
    ## LCM_3  1.0000000 0.6473983 0.5666346
    ## LCM_32 0.6473983 1.0000000 0.6036661
    ## LCM_33 0.5666346 0.6036661 1.0000000

    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter
    ## Warning in par(usr): argument 1 does not name a graphical parameter

    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'
    ## Warning in KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, :
    ## Binning grid too coarse for current (small) bandwidth: consider increasing
    ## 'gridsize'

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### 0.3.2 batch effects

``` r
#as=assocComp(mBase=meth_filter_destrand,dplyr::select(meta,c("PCR_ReAmp_Cycles", "Fragment")))
#as
```

### 0.3.3 other possible filtering

``` r
# get percent methylation matrix
pm=percMethylation(meth_filter_destrand)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)
```

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# keep only CpG with standard deviations larger than 2%
#meth <- meth_filter_destrand[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth_filter_destrand)
```

    ## [1] 134034

``` r
#nrow(meth)
```

### 0.3.4 Identify DML

``` r
DMLStats_Tissue <- methylKit::calculateDiffMeth(meth_filter_destrand, overdispersion = "MN", test = "Chisq", mc.cores = 8) #Calculate differential methylation statistics and include covariate information.
```

    ## two groups detected:
    ##  will calculate methylation difference as the difference of
    ## treatment (group: 1) - control (group: 0)

``` r
head(DMLStats_Tissue) #Look at differential methylation output
```

    ##                                  chr  start    end strand      pvalue    qvalue
    ## 1 Pocillopora_acuta_HIv2___Sc0000000   5326   5326      + 1.000000000 0.8599801
    ## 2 Pocillopora_acuta_HIv2___Sc0000000   5331   5331      + 1.000000000 0.8599801
    ## 3 Pocillopora_acuta_HIv2___Sc0000000  83557  83557      + 0.004626153 0.5025833
    ## 4 Pocillopora_acuta_HIv2___Sc0000000  83565  83565      + 0.053009864 0.8599801
    ## 5 Pocillopora_acuta_HIv2___Sc0000000  83588  83588      + 0.123774009 0.8599801
    ## 6 Pocillopora_acuta_HIv2___Sc0000000 153088 153088      + 1.000000000 0.8599801
    ##   meth.diff
    ## 1   0.00000
    ## 2   0.00000
    ## 3 -54.54545
    ## 4 -28.12500
    ## 5 -26.66667
    ## 6   0.00000

``` r
# Filter DMRs with q-value < 0.05
significant_dmg <- getData(DMLStats_Tissue[DMLStats_Tissue$qvalue < 0.05, ])

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

``` r
DMLs <- methylKit::getMethylDiff(DMLStats_Tissue, difference = 2, qvalue = 0.05) #Identify DML based on difference threshold

length(DMLs$chr) #DML
```

    ## [1] 190

``` r
head(DMLs)
```

    ##                                     chr   start     end strand       pvalue
    ## 4981 Pocillopora_acuta_HIv2___Sc0000002 3933210 3933210      + 6.801995e-06
    ## 5079 Pocillopora_acuta_HIv2___Sc0000002 4354140 4354140      + 4.217920e-05
    ## 5080 Pocillopora_acuta_HIv2___Sc0000002 4354519 4354519      + 3.498027e-06
    ## 5083 Pocillopora_acuta_HIv2___Sc0000002 4354790 4354790      + 2.793794e-10
    ## 6087 Pocillopora_acuta_HIv2___Sc0000002 9281328 9281328      + 5.825224e-05
    ## 7904 Pocillopora_acuta_HIv2___Sc0000003 2149671 2149671      + 1.223095e-06
    ##            qvalue meth.diff
    ## 4981 8.340879e-03  75.00000
    ## 5079 3.077121e-02  15.29412
    ## 5080 5.236437e-03  22.64151
    ## 5083 2.477162e-06  29.78723
    ## 6087 3.837845e-02  23.07692
    ## 7904 2.468497e-03 -44.00000

## 0.4 Further look at genome wide methylation

``` r
diffMethPerChr(DMLStats_Tissue, meth.cutoff = 2, qvalue.cutoff = 0.05,cex.names=.75)
```

    ## Warning in eval(quote(list(...)), env): NAs introduced by coercion

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

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
DML_grange
```

    ## GRanges object with 190 ranges and 3 metadata columns:
    ##                       seqnames    ranges strand |      pvalue      qvalue
    ##                          <Rle> <IRanges>  <Rle> |   <numeric>   <numeric>
    ##     [1] Pocillopora_acuta_HI..   3933210      + | 6.80199e-06 8.34088e-03
    ##     [2] Pocillopora_acuta_HI..   4354140      + | 4.21792e-05 3.07712e-02
    ##     [3] Pocillopora_acuta_HI..   4354519      + | 3.49803e-06 5.23644e-03
    ##     [4] Pocillopora_acuta_HI..   4354790      + | 2.79379e-10 2.47716e-06
    ##     [5] Pocillopora_acuta_HI..   9281328      + | 5.82522e-05 3.83785e-02
    ##     ...                    ...       ...    ... .         ...         ...
    ##   [186] Pocillopora_acuta_HI..    442312      + | 3.87231e-07 0.001062733
    ##   [187] Pocillopora_acuta_HI..    442314      + | 4.90440e-07 0.001256253
    ##   [188] Pocillopora_acuta_HI..    893636      + | 3.35911e-08 0.000147786
    ##   [189] Pocillopora_acuta_HI..    903270      + | 7.10337e-05 0.043785082
    ##   [190] Pocillopora_acuta_HI..    903292      + | 7.07557e-05 0.043785082
    ##         meth.diff
    ##         <numeric>
    ##     [1]   75.0000
    ##     [2]   15.2941
    ##     [3]   22.6415
    ##     [4]   29.7872
    ##     [5]   23.0769
    ##     ...       ...
    ##   [186]   47.9526
    ##   [187]   47.8618
    ##   [188]   69.4872
    ##   [189]   44.4444
    ##   [190]   40.8192
    ##   -------
    ##   seqinfo: 64 sequences from an unspecified genome; no seqlengths

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

    ##                              DML_chr DML_start DML_end   DML_qvalue
    ## 1 Pocillopora_acuta_HIv2___Sc0000002   4354140 4354140 3.077121e-02
    ## 2 Pocillopora_acuta_HIv2___Sc0000002   4354519 4354519 5.236437e-03
    ## 3 Pocillopora_acuta_HIv2___Sc0000002   4354790 4354790 2.477162e-06
    ## 4 Pocillopora_acuta_HIv2___Sc0000002   9281328 9281328 3.837845e-02
    ## 5 Pocillopora_acuta_HIv2___Sc0000003   2461240 2461240 3.643553e-02
    ## 6 Pocillopora_acuta_HIv2___Sc0000003   8702392 8702392 1.725322e-02
    ##   DML_methdiff                     transcript_chr transcript_start
    ## 1     15.29412 Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 2     22.64151 Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 3     29.78723 Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 4     23.07692 Pocillopora_acuta_HIv2___Sc0000002          9279157
    ## 5    -35.00000 Pocillopora_acuta_HIv2___Sc0000003          2441368
    ## 6     30.30303 Pocillopora_acuta_HIv2___Sc0000003          8696824
    ##   transcript_end                             transcript_id
    ## 1        4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 2        4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 3        4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 4        9293249 Pocillopora_acuta_HIv2___RNAseq.g25453.t1
    ## 5        2461336  Pocillopora_acuta_HIv2___RNAseq.g4860.t2
    ## 6        8734794  Pocillopora_acuta_HIv2___RNAseq.g5465.t1
    ##                                     gene_id
    ## 1 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 2 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 3 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 4 Pocillopora_acuta_HIv2___RNAseq.g25453.t1
    ## 5  Pocillopora_acuta_HIv2___RNAseq.g4860.t2
    ## 6  Pocillopora_acuta_HIv2___RNAseq.g5465.t1

``` r
# how many DMLs are in gene bodies?
length(DML_transcript_annot$gene_id)
```

    ## [1] 91

``` r
#how many genes do these consist of?
length(unique(DML_transcript_annot$gene_id))
```

    ## [1] 65

``` r
DML_transcript_annot %>% group_by(transcript_id) %>% summarize(num_DMLs=n()) %>% summary()
```

    ##  transcript_id         num_DMLs  
    ##  Length:65          Min.   :1.0  
    ##  Class :character   1st Qu.:1.0  
    ##  Mode  :character   Median :1.0  
    ##                     Mean   :1.4  
    ##                     3rd Qu.:1.0  
    ##                     Max.   :7.0

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

    ##                                 DML_chr DML_start DML_end   DML_qvalue
    ## 1    Pocillopora_acuta_HIv2___Sc0000002   4354140 4354140 3.077121e-02
    ## 2    Pocillopora_acuta_HIv2___Sc0000002   4354519 4354519 5.236437e-03
    ## 3    Pocillopora_acuta_HIv2___Sc0000002   4354790 4354790 2.477162e-06
    ## 10   Pocillopora_acuta_HIv2___Sc0000006   5260722 5260722 2.225195e-02
    ## 11   Pocillopora_acuta_HIv2___Sc0000007   1469599 1469599 1.443216e-02
    ## 19   Pocillopora_acuta_HIv2___Sc0000009   1926284 1926284 4.300753e-02
    ## 20   Pocillopora_acuta_HIv2___Sc0000009   1926297 1926297 3.564570e-02
    ## 37 Pocillopora_acuta_HIv2___xfSc0000000   5400485 5400485 3.281083e-02
    ## 63 Pocillopora_acuta_HIv2___xfSc0000009   1377935 1377935 1.383141e-02
    ## 64 Pocillopora_acuta_HIv2___xfSc0000009   3043722 3043722 1.327307e-02
    ## 74 Pocillopora_acuta_HIv2___xfSc0000024   2471817 2471817 2.468497e-03
    ##    DML_methdiff                       transcript_chr transcript_start
    ## 1     15.294118   Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 2     22.641509   Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 3     29.787234   Pocillopora_acuta_HIv2___Sc0000002          4352529
    ## 10    36.502508   Pocillopora_acuta_HIv2___Sc0000006          5258672
    ## 11   -29.268293   Pocillopora_acuta_HIv2___Sc0000007          1469481
    ## 19    15.789474   Pocillopora_acuta_HIv2___Sc0000009          1924903
    ## 20    14.754098   Pocillopora_acuta_HIv2___Sc0000009          1924903
    ## 37    25.806452 Pocillopora_acuta_HIv2___xfSc0000000          5393991
    ## 63    14.035088 Pocillopora_acuta_HIv2___xfSc0000009          1376679
    ## 64   -55.555556 Pocillopora_acuta_HIv2___xfSc0000009          3018389
    ## 74    -9.360961 Pocillopora_acuta_HIv2___xfSc0000024          2471129
    ##    transcript_end                             transcript_id
    ## 1         4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 2         4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 3         4355343 Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 10        5262000 Pocillopora_acuta_HIv2___RNAseq.g15474.t1
    ## 11        1470395 Pocillopora_acuta_HIv2___RNAseq.g17564.t1
    ## 19        1932130  Pocillopora_acuta_HIv2___RNAseq.g2626.t1
    ## 20        1932130  Pocillopora_acuta_HIv2___RNAseq.g2626.t1
    ## 37        5402276 Pocillopora_acuta_HIv2___RNAseq.g22998.t1
    ## 63        1378770 Pocillopora_acuta_HIv2___RNAseq.g12189.t1
    ## 64        3057287      Pocillopora_acuta_HIv2___TS.g3256.t1
    ## 74        2473563 Pocillopora_acuta_HIv2___RNAseq.g20860.t1
    ##                                      gene_id
    ## 1  Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 2  Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 3  Pocillopora_acuta_HIv2___RNAseq.g25038.t1
    ## 10 Pocillopora_acuta_HIv2___RNAseq.g15474.t1
    ## 11 Pocillopora_acuta_HIv2___RNAseq.g17564.t1
    ## 19  Pocillopora_acuta_HIv2___RNAseq.g2626.t1
    ## 20  Pocillopora_acuta_HIv2___RNAseq.g2626.t1
    ## 37 Pocillopora_acuta_HIv2___RNAseq.g22998.t1
    ## 63 Pocillopora_acuta_HIv2___RNAseq.g12189.t1
    ## 64      Pocillopora_acuta_HIv2___TS.g3256.t1
    ## 74 Pocillopora_acuta_HIv2___RNAseq.g20860.t1

``` r
DE_05[DE_05$query %in% DML_transcript_annot$transcript_id,]
```

    ##                                          query   baseMean log2FoldChange
    ## 21   Pocillopora_acuta_HIv2___RNAseq.g20860.t1 1265.27422    -15.5762448
    ## 1788 Pocillopora_acuta_HIv2___RNAseq.g17564.t1  708.27278     -0.7085281
    ## 1802  Pocillopora_acuta_HIv2___RNAseq.g2626.t1  237.04137      5.6905130
    ## 1928 Pocillopora_acuta_HIv2___RNAseq.g15474.t1   42.82082     -4.5903434
    ## 2005 Pocillopora_acuta_HIv2___RNAseq.g22998.t1  287.46168     -0.7535547
    ## 2924 Pocillopora_acuta_HIv2___RNAseq.g25038.t1  453.44869      1.4528771
    ## 3133 Pocillopora_acuta_HIv2___RNAseq.g12189.t1  184.60827     -0.3170635
    ## 3249      Pocillopora_acuta_HIv2___TS.g3256.t1  295.29556     -1.0719805
    ##         lfcSE       pvalue         padj
    ## 21   3.873810 1.358173e-16 9.354577e-14
    ## 1788 1.285995 1.192290e-04 9.645014e-04
    ## 1802 2.267624 1.246384e-04 1.000427e-03
    ## 1928 2.195922 1.899619e-04 1.425108e-03
    ## 2005 1.430480 2.545588e-04 1.836378e-03
    ## 2924 2.114543 3.385288e-03 1.674583e-02
    ## 3133 1.047937 5.237407e-03 2.417933e-02
    ## 3249 1.600315 6.755427e-03 3.007402e-02

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

    ## Warning: ggrepel: 3 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# how many CpG sites are there after filtering?
dim(meth_filter_destrand)
```

    ## [1] 134034     34

``` r
ML_grange = as(meth_filter_destrand,"GRanges")

# identify all the MLs that are in transcripts
genes_with_MLs <- findOverlaps(ML_grange, transcripts,ignore.strand = TRUE)
```

    ## Warning in .merge_two_Seqinfo_objects(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': Pocillopora_acuta_HIv2___Sc0000054, Pocillopora_acuta_HIv2___Sc0000069, Pocillopora_acuta_HIv2___xfSc0000067, Pocillopora_acuta_HIv2___xfSc0000112, Pocillopora_acuta_HIv2___xfSc0000153, Pocillopora_acuta_HIv2___xfSc0000172, Pocillopora_acuta_HIv2___xfSc0000225, Pocillopora_acuta_HIv2___xfSc0000365, Pocillopora_acuta_HIv2___xfSc0000369, Pocillopora_acuta_HIv2___xfSc0000372, Pocillopora_acuta_HIv2___xfSc0000375, Pocillopora_acuta_HIv2___xfSc0000383, Pocillopora_acuta_HIv2___xpSc0000428
    ##   - in 'y': Pocillopora_acuta_HIv2___Sc0000056, Pocillopora_acuta_HIv2___Sc0000058, Pocillopora_acuta_HIv2___Sc0000059, Pocillopora_acuta_HIv2___Sc0000062, Pocillopora_acuta_HIv2___Sc0000064, Pocillopora_acuta_HIv2___Sc0000068, Pocillopora_acuta_HIv2___Sc0000070, Pocillopora_acuta_HIv2___Sc0000073, Pocillopora_acuta_HIv2___Sc0000076, Pocillopora_acuta_HIv2___Sc0000077, Pocillopora_acuta_HIv2___Sc0000080, Pocillopora_acuta_HIv2___xfSc0000071, Pocillopora_acuta_HIv2___xfSc0000074, Pocillopora_acuta_HIv2___xfSc0000078, Pocillopora_acuta_HIv2___xfSc0000087, Pocillopora_acuta_HIv2___xfSc0000089, Pocillopora_acuta_HIv2___xfSc0000090, Pocillopora_acuta_HIv2___xfSc0000095, Pocillopora_acuta_HIv2___xfSc0000097, Pocillopora_acuta_HIv2___xfSc0000099, Pocillopora_acuta_HIv2___xfSc0000104, Pocillopora_acuta_HIv2___xfSc0000105, Pocillopora_acuta_HIv2___xfSc0000108, Pocillopora_acuta_HIv2___xfSc0000109, Pocillopora_acuta_HIv2___xfSc0000110, Pocillopora_acuta_HIv2___xfSc0000116, Pocillopora_acuta_HIv2___xfSc0000119, Pocillopora_acuta_HIv2___xfSc0000121, Pocillopora_acuta_HIv2___xfSc0000122, Pocillopora_acuta_HIv2___xfSc0000128, Pocillopora_acuta_HIv2___xfSc0000129, Pocillopora_acuta_HIv2___xfSc0000131, Pocillopora_acuta_HIv2___xfSc0000132, Pocillopora_acuta_HIv2___xfSc0000135, Pocillopora_acuta_HIv2___xfSc0000136, Pocillopora_acuta_HIv2___xfSc0000139, Pocillopora_acuta_HIv2___xfSc0000140, Pocillopora_acuta_HIv2___xfSc0000141, Pocillopora_acuta_HIv2___xfSc0000143, Pocillopora_acuta_HIv2___xfSc0000144, Pocillopora_acuta_HIv2___xfSc0000147, Pocillopora_acuta_HIv2___xfSc0000148, Pocillopora_acuta_HIv2___xfSc0000151, Pocillopora_acuta_HIv2___xfSc0000152, Pocillopora_acuta_HIv2___xfSc0000154, Pocillopora_acuta_HIv2___xfSc0000155, Pocillopora_acuta_HIv2___xfSc0000158, Pocillopora_acuta_HIv2___xfSc0000160, Pocillopora_acuta_HIv2___xfSc0000161, Pocillopora_acuta_HIv2___xfSc0000163, Pocillopora_acuta_HIv2___xfSc0000165, Pocillopora_acuta_HIv2___xfSc0000166, Pocillopora_acuta_HIv2___xfSc0000167, Pocillopora_acuta_HIv2___xfSc0000168, Pocillopora_acuta_HIv2___xfSc0000171, Pocillopora_acuta_HIv2___xfSc0000177, Pocillopora_acuta_HIv2___xfSc0000179, Pocillopora_acuta_HIv2___xfSc0000180, Pocillopora_acuta_HIv2___xfSc0000182, Pocillopora_acuta_HIv2___xfSc0000183, Pocillopora_acuta_HIv2___xfSc0000184, Pocillopora_acuta_HIv2___xfSc0000185, Pocillopora_acuta_HIv2___xfSc0000186, Pocillopora_acuta_HIv2___xfSc0000187, Pocillopora_acuta_HIv2___xfSc0000189, Pocillopora_acuta_HIv2___xfSc0000191, Pocillopora_acuta_HIv2___xfSc0000192, Pocillopora_acuta_HIv2___xfSc0000193, Pocillopora_acuta_HIv2___xfSc0000195, Pocillopora_acuta_HIv2___xfSc0000196, Pocillopora_acuta_HIv2___xfSc0000200, Pocillopora_acuta_HIv2___xfSc0000203, Pocillopora_acuta_HIv2___xfSc0000204, Pocillopora_acuta_HIv2___xfSc0000205, Pocillopora_acuta_HIv2___xfSc0000207, Pocillopora_acuta_HIv2___xfSc0000208, Pocillopora_acuta_HIv2___xfSc0000209, Pocillopora_acuta_HIv2___xfSc0000212, Pocillopora_acuta_HIv2___xfSc0000214, Pocillopora_acuta_HIv2___xfSc0000215, Pocillopora_acuta_HIv2___xfSc0000216, Pocillopora_acuta_HIv2___xfSc0000217, Pocillopora_acuta_HIv2___xfSc0000218, Pocillopora_acuta_HIv2___xfSc0000223, Pocillopora_acuta_HIv2___xfSc0000228, Pocillopora_acuta_HIv2___xfSc0000230, Pocillopora_acuta_HIv2___xfSc0000231, Pocillopora_acuta_HIv2___xfSc0000232, Pocillopora_acuta_HIv2___xfSc0000233, Pocillopora_acuta_HIv2___xfSc0000234, Pocillopora_acuta_HIv2___xfSc0000236, Pocillopora_acuta_HIv2___xfSc0000237, Pocillopora_acuta_HIv2___xfSc0000239, Pocillopora_acuta_HIv2___xfSc0000240, Pocillopora_acuta_HIv2___xfSc0000243, Pocillopora_acuta_HIv2___xfSc0000244, Pocillopora_acuta_HIv2___xfSc0000248, Pocillopora_acuta_HIv2___xfSc0000249, Pocillopora_acuta_HIv2___xfSc0000250, Pocillopora_acuta_HIv2___xfSc0000251, Pocillopora_acuta_HIv2___xfSc0000252, Pocillopora_acuta_HIv2___xfSc0000255, Pocillopora_acuta_HIv2___xfSc0000256, Pocillopora_acuta_HIv2___xfSc0000257, Pocillopora_acuta_HIv2___xfSc0000258, Pocillopora_acuta_HIv2___xfSc0000260, Pocillopora_acuta_HIv2___xfSc0000262, Pocillopora_acuta_HIv2___xfSc0000264, Pocillopora_acuta_HIv2___xfSc0000265, Pocillopora_acuta_HIv2___xfSc0000266, Pocillopora_acuta_HIv2___xfSc0000267, Pocillopora_acuta_HIv2___xfSc0000268, Pocillopora_acuta_HIv2___xfSc0000269, Pocillopora_acuta_HIv2___xfSc0000271, Pocillopora_acuta_HIv2___xfSc0000272, Pocillopora_acuta_HIv2___xfSc0000273, Pocillopora_acuta_HIv2___xfSc0000275, Pocillopora_acuta_HIv2___xfSc0000276, Pocillopora_acuta_HIv2___xfSc0000277, Pocillopora_acuta_HIv2___xfSc0000280, Pocillopora_acuta_HIv2___xfSc0000281, Pocillopora_acuta_HIv2___xfSc0000283, Pocillopora_acuta_HIv2___xfSc0000284, Pocillopora_acuta_HIv2___xfSc0000287, Pocillopora_acuta_HIv2___xfSc0000288, Pocillopora_acuta_HIv2___xfSc0000291, Pocillopora_acuta_HIv2___xfSc0000293, Pocillopora_acuta_HIv2___xfSc0000296, Pocillopora_acuta_HIv2___xfSc0000302, Pocillopora_acuta_HIv2___xfSc0000303, Pocillopora_acuta_HIv2___xfSc0000305, Pocillopora_acuta_HIv2___xfSc0000306, Pocillopora_acuta_HIv2___xfSc0000310, Pocillopora_acuta_HIv2___xfSc0000311, Pocillopora_acuta_HIv2___xfSc0000313, Pocillopora_acuta_HIv2___xfSc0000315, Pocillopora_acuta_HIv2___xfSc0000316, Pocillopora_acuta_HIv2___xfSc0000319, Pocillopora_acuta_HIv2___xfSc0000322, Pocillopora_acuta_HIv2___xfSc0000325, Pocillopora_acuta_HIv2___xfSc0000326, Pocillopora_acuta_HIv2___xfSc0000327, Pocillopora_acuta_HIv2___xfSc0000334, Pocillopora_acuta_HIv2___xfSc0000343, Pocillopora_acuta_HIv2___xfSc0000344, Pocillopora_acuta_HIv2___xfSc0000345, Pocillopora_acuta_HIv2___xfSc0000346, Pocillopora_acuta_HIv2___xfSc0000347, Pocillopora_acuta_HIv2___xfSc0000349, Pocillopora_acuta_HIv2___xfSc0000350, Pocillopora_acuta_HIv2___xfSc0000354, Pocillopora_acuta_HIv2___xfSc0000361, Pocillopora_acuta_HIv2___xfSc0000362, Pocillopora_acuta_HIv2___xfSc0000363, Pocillopora_acuta_HIv2___xfSc0000368, Pocillopora_acuta_HIv2___xfSc0000370, Pocillopora_acuta_HIv2___xfSc0000371, Pocillopora_acuta_HIv2___xfSc0000374, Pocillopora_acuta_HIv2___xfSc0000376, Pocillopora_acuta_HIv2___xfSc0000377, Pocillopora_acuta_HIv2___xfSc0000378, Pocillopora_acuta_HIv2___xfSc0000379, Pocillopora_acuta_HIv2___xfSc0000380, Pocillopora_acuta_HIv2___xfSc0000384, Pocillopora_acuta_HIv2___xfSc0000386, Pocillopora_acuta_HIv2___xpSc0000402, Pocillopora_acuta_HIv2___xpSc0000412, Pocillopora_acuta_HIv2___xpSc0000414, Pocillopora_acuta_HIv2___xpSc0000415, Pocillopora_acuta_HIv2___xpSc0000417, Pocillopora_acuta_HIv2___xpSc0000422, Pocillopora_acuta_HIv2___xpSc0000423, Pocillopora_acuta_HIv2___xpSc0000425, Pocillopora_acuta_HIv2___xpSc0000430, Pocillopora_acuta_HIv2___xpSc0000439, Pocillopora_acuta_HIv2___xpSc0000447
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

``` r
# Extract matching transcript information
ML_transcript_annot <- data.frame(
  ML_chr = seqnames(ML_grange)[queryHits(genes_with_MLs)],
  ML_start = start(ML_grange)[queryHits(genes_with_MLs)],
  ML_end = end(ML_grange)[queryHits(genes_with_MLs)],
  transcript_chr = seqnames(transcripts)[subjectHits(genes_with_MLs)],
  transcript_start = start(transcripts)[subjectHits(genes_with_MLs)],
  transcript_end = end(transcripts)[subjectHits(genes_with_MLs)],
  transcript_id = transcripts$transcript_id[subjectHits(genes_with_MLs)],
  gene_id = transcripts$gene_id[subjectHits(genes_with_MLs)]
)

#how many MLs are represented
length(ML_transcript_annot$transcript_id)
```

    ## [1] 38747

``` r
# Find overlaps between methylated loci and transcripts, averaged by gene
MLs_in_genes <- regionCounts(meth_filter_destrand, regions=transcripts)
```

    ## Warning in .merge_two_Seqinfo_objects(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': Pocillopora_acuta_HIv2___Sc0000056, Pocillopora_acuta_HIv2___Sc0000058, Pocillopora_acuta_HIv2___Sc0000059, Pocillopora_acuta_HIv2___Sc0000062, Pocillopora_acuta_HIv2___Sc0000064, Pocillopora_acuta_HIv2___Sc0000068, Pocillopora_acuta_HIv2___Sc0000070, Pocillopora_acuta_HIv2___Sc0000073, Pocillopora_acuta_HIv2___Sc0000076, Pocillopora_acuta_HIv2___Sc0000077, Pocillopora_acuta_HIv2___Sc0000080, Pocillopora_acuta_HIv2___xfSc0000071, Pocillopora_acuta_HIv2___xfSc0000074, Pocillopora_acuta_HIv2___xfSc0000078, Pocillopora_acuta_HIv2___xfSc0000087, Pocillopora_acuta_HIv2___xfSc0000089, Pocillopora_acuta_HIv2___xfSc0000090, Pocillopora_acuta_HIv2___xfSc0000095, Pocillopora_acuta_HIv2___xfSc0000097, Pocillopora_acuta_HIv2___xfSc0000099, Pocillopora_acuta_HIv2___xfSc0000104, Pocillopora_acuta_HIv2___xfSc0000105, Pocillopora_acuta_HIv2___xfSc0000108, Pocillopora_acuta_HIv2___xfSc0000109, Pocillopora_acuta_HIv2___xfSc0000110, Pocillopora_acuta_HIv2___xfSc0000116, Pocillopora_acuta_HIv2___xfSc0000119, Pocillopora_acuta_HIv2___xfSc0000121, Pocillopora_acuta_HIv2___xfSc0000122, Pocillopora_acuta_HIv2___xfSc0000128, Pocillopora_acuta_HIv2___xfSc0000129, Pocillopora_acuta_HIv2___xfSc0000131, Pocillopora_acuta_HIv2___xfSc0000132, Pocillopora_acuta_HIv2___xfSc0000135, Pocillopora_acuta_HIv2___xfSc0000136, Pocillopora_acuta_HIv2___xfSc0000139, Pocillopora_acuta_HIv2___xfSc0000140, Pocillopora_acuta_HIv2___xfSc0000141, Pocillopora_acuta_HIv2___xfSc0000143, Pocillopora_acuta_HIv2___xfSc0000144, Pocillopora_acuta_HIv2___xfSc0000147, Pocillopora_acuta_HIv2___xfSc0000148, Pocillopora_acuta_HIv2___xfSc0000151, Pocillopora_acuta_HIv2___xfSc0000152, Pocillopora_acuta_HIv2___xfSc0000154, Pocillopora_acuta_HIv2___xfSc0000155, Pocillopora_acuta_HIv2___xfSc0000158, Pocillopora_acuta_HIv2___xfSc0000160, Pocillopora_acuta_HIv2___xfSc0000161, Pocillopora_acuta_HIv2___xfSc0000163, Pocillopora_acuta_HIv2___xfSc0000165, Pocillopora_acuta_HIv2___xfSc0000166, Pocillopora_acuta_HIv2___xfSc0000167, Pocillopora_acuta_HIv2___xfSc0000168, Pocillopora_acuta_HIv2___xfSc0000171, Pocillopora_acuta_HIv2___xfSc0000177, Pocillopora_acuta_HIv2___xfSc0000179, Pocillopora_acuta_HIv2___xfSc0000180, Pocillopora_acuta_HIv2___xfSc0000182, Pocillopora_acuta_HIv2___xfSc0000183, Pocillopora_acuta_HIv2___xfSc0000184, Pocillopora_acuta_HIv2___xfSc0000185, Pocillopora_acuta_HIv2___xfSc0000186, Pocillopora_acuta_HIv2___xfSc0000187, Pocillopora_acuta_HIv2___xfSc0000189, Pocillopora_acuta_HIv2___xfSc0000191, Pocillopora_acuta_HIv2___xfSc0000192, Pocillopora_acuta_HIv2___xfSc0000193, Pocillopora_acuta_HIv2___xfSc0000195, Pocillopora_acuta_HIv2___xfSc0000196, Pocillopora_acuta_HIv2___xfSc0000200, Pocillopora_acuta_HIv2___xfSc0000203, Pocillopora_acuta_HIv2___xfSc0000204, Pocillopora_acuta_HIv2___xfSc0000205, Pocillopora_acuta_HIv2___xfSc0000207, Pocillopora_acuta_HIv2___xfSc0000208, Pocillopora_acuta_HIv2___xfSc0000209, Pocillopora_acuta_HIv2___xfSc0000212, Pocillopora_acuta_HIv2___xfSc0000214, Pocillopora_acuta_HIv2___xfSc0000215, Pocillopora_acuta_HIv2___xfSc0000216, Pocillopora_acuta_HIv2___xfSc0000217, Pocillopora_acuta_HIv2___xfSc0000218, Pocillopora_acuta_HIv2___xfSc0000223, Pocillopora_acuta_HIv2___xfSc0000228, Pocillopora_acuta_HIv2___xfSc0000230, Pocillopora_acuta_HIv2___xfSc0000231, Pocillopora_acuta_HIv2___xfSc0000232, Pocillopora_acuta_HIv2___xfSc0000233, Pocillopora_acuta_HIv2___xfSc0000234, Pocillopora_acuta_HIv2___xfSc0000236, Pocillopora_acuta_HIv2___xfSc0000237, Pocillopora_acuta_HIv2___xfSc0000239, Pocillopora_acuta_HIv2___xfSc0000240, Pocillopora_acuta_HIv2___xfSc0000243, Pocillopora_acuta_HIv2___xfSc0000244, Pocillopora_acuta_HIv2___xfSc0000248, Pocillopora_acuta_HIv2___xfSc0000249, Pocillopora_acuta_HIv2___xfSc0000250, Pocillopora_acuta_HIv2___xfSc0000251, Pocillopora_acuta_HIv2___xfSc0000252, Pocillopora_acuta_HIv2___xfSc0000255, Pocillopora_acuta_HIv2___xfSc0000256, Pocillopora_acuta_HIv2___xfSc0000257, Pocillopora_acuta_HIv2___xfSc0000258, Pocillopora_acuta_HIv2___xfSc0000260, Pocillopora_acuta_HIv2___xfSc0000262, Pocillopora_acuta_HIv2___xfSc0000264, Pocillopora_acuta_HIv2___xfSc0000265, Pocillopora_acuta_HIv2___xfSc0000266, Pocillopora_acuta_HIv2___xfSc0000267, Pocillopora_acuta_HIv2___xfSc0000268, Pocillopora_acuta_HIv2___xfSc0000269, Pocillopora_acuta_HIv2___xfSc0000271, Pocillopora_acuta_HIv2___xfSc0000272, Pocillopora_acuta_HIv2___xfSc0000273, Pocillopora_acuta_HIv2___xfSc0000275, Pocillopora_acuta_HIv2___xfSc0000276, Pocillopora_acuta_HIv2___xfSc0000277, Pocillopora_acuta_HIv2___xfSc0000280, Pocillopora_acuta_HIv2___xfSc0000281, Pocillopora_acuta_HIv2___xfSc0000283, Pocillopora_acuta_HIv2___xfSc0000284, Pocillopora_acuta_HIv2___xfSc0000287, Pocillopora_acuta_HIv2___xfSc0000288, Pocillopora_acuta_HIv2___xfSc0000291, Pocillopora_acuta_HIv2___xfSc0000293, Pocillopora_acuta_HIv2___xfSc0000296, Pocillopora_acuta_HIv2___xfSc0000302, Pocillopora_acuta_HIv2___xfSc0000303, Pocillopora_acuta_HIv2___xfSc0000305, Pocillopora_acuta_HIv2___xfSc0000306, Pocillopora_acuta_HIv2___xfSc0000310, Pocillopora_acuta_HIv2___xfSc0000311, Pocillopora_acuta_HIv2___xfSc0000313, Pocillopora_acuta_HIv2___xfSc0000315, Pocillopora_acuta_HIv2___xfSc0000316, Pocillopora_acuta_HIv2___xfSc0000319, Pocillopora_acuta_HIv2___xfSc0000322, Pocillopora_acuta_HIv2___xfSc0000325, Pocillopora_acuta_HIv2___xfSc0000326, Pocillopora_acuta_HIv2___xfSc0000327, Pocillopora_acuta_HIv2___xfSc0000334, Pocillopora_acuta_HIv2___xfSc0000343, Pocillopora_acuta_HIv2___xfSc0000344, Pocillopora_acuta_HIv2___xfSc0000345, Pocillopora_acuta_HIv2___xfSc0000346, Pocillopora_acuta_HIv2___xfSc0000347, Pocillopora_acuta_HIv2___xfSc0000349, Pocillopora_acuta_HIv2___xfSc0000350, Pocillopora_acuta_HIv2___xfSc0000354, Pocillopora_acuta_HIv2___xfSc0000361, Pocillopora_acuta_HIv2___xfSc0000362, Pocillopora_acuta_HIv2___xfSc0000363, Pocillopora_acuta_HIv2___xfSc0000368, Pocillopora_acuta_HIv2___xfSc0000370, Pocillopora_acuta_HIv2___xfSc0000371, Pocillopora_acuta_HIv2___xfSc0000374, Pocillopora_acuta_HIv2___xfSc0000376, Pocillopora_acuta_HIv2___xfSc0000377, Pocillopora_acuta_HIv2___xfSc0000378, Pocillopora_acuta_HIv2___xfSc0000379, Pocillopora_acuta_HIv2___xfSc0000380, Pocillopora_acuta_HIv2___xfSc0000384, Pocillopora_acuta_HIv2___xfSc0000386, Pocillopora_acuta_HIv2___xpSc0000402, Pocillopora_acuta_HIv2___xpSc0000412, Pocillopora_acuta_HIv2___xpSc0000414, Pocillopora_acuta_HIv2___xpSc0000415, Pocillopora_acuta_HIv2___xpSc0000417, Pocillopora_acuta_HIv2___xpSc0000422, Pocillopora_acuta_HIv2___xpSc0000423, Pocillopora_acuta_HIv2___xpSc0000425, Pocillopora_acuta_HIv2___xpSc0000430, Pocillopora_acuta_HIv2___xpSc0000439, Pocillopora_acuta_HIv2___xpSc0000447
    ##   - in 'y': Pocillopora_acuta_HIv2___Sc0000054, Pocillopora_acuta_HIv2___Sc0000069, Pocillopora_acuta_HIv2___xfSc0000067, Pocillopora_acuta_HIv2___xfSc0000112, Pocillopora_acuta_HIv2___xfSc0000153, Pocillopora_acuta_HIv2___xfSc0000172, Pocillopora_acuta_HIv2___xfSc0000225, Pocillopora_acuta_HIv2___xfSc0000365, Pocillopora_acuta_HIv2___xfSc0000369, Pocillopora_acuta_HIv2___xfSc0000372, Pocillopora_acuta_HIv2___xfSc0000375, Pocillopora_acuta_HIv2___xfSc0000383, Pocillopora_acuta_HIv2___xpSc0000428
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

``` r
MLs_in_genes$gene_id <- transcripts$gene_id[match(paste(MLs_in_genes$chr, MLs_in_genes$start, MLs_in_genes$end), paste(seqnames(transcripts), start(transcripts), end(transcripts)))]

#how many genes are represented
nrow(MLs_in_genes)
```

    ## [1] 4972

``` r
#does this match the above?
length(unique(ML_transcript_annot$transcript_id))
```

    ## [1] 4972

``` r
percent_meth <- percMethylation(MLs_in_genes)
percent_meth <- as.data.frame(percent_meth)
percent_meth$gene_id <- MLs_in_genes$gene_id 
percent_meth <- percent_meth %>% select(gene_id, everything())

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

Now, combine the genes that contain methylated loci that passed
filtering with the genes that were found to be expressed in our RNA-seq
dataset (14464).

``` r
# number of genes containing MLs
nrow(percent_meth)
```

    ## [1] 4972

``` r
# number of expressed genes
nrow(DESeq)
```

    ## [1] 14464

``` r
expressed_genes_percent_meth <- merge(percent_meth, DESeq, by.x = "gene_id", by.y = "query")

nrow(expressed_genes_percent_meth)
```

    ## [1] 2572

``` r
plot_data <- merge(percent_meth_long, DESeq, by.x = "gene_id", by.y = "query")
```

### 0.5.1 BaseMean: Raw gene expression levels (skewed by highly expressed genes)

``` r
ggplot(plot_data, aes(y = baseMean, x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "DeSeq2 BaseMean", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
ggplot(plot_data, aes(y = baseMean, x = tissue_percent_meth,color=Tissue)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "DeSeq2 BaseMean", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

### 0.5.2 Log2 of BaseMean: Transformed gene expression levels

``` r
ggplot(plot_data, aes(y = log2(baseMean), x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Log2(DeSeq2 BaseMean)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
ggplot(plot_data, aes(y = log2(baseMean), x = tissue_percent_meth,color=Tissue)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Log2(DeSeq2 BaseMean)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

### 0.5.3 Absolute value of Log2FoldChange: Relationship between overall methylation of a gene and whether or not it is differentially expressed

``` r
ggplot(plot_data, aes(y = abs(log2FoldChange), x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +  stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Abs(DeSeq2 Log2FoldChange)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
ggplot(plot_data, aes(y = abs(log2FoldChange), x = tissue_percent_meth,color=Tissue)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +  stat_poly_eq(use_label("eq", "R2")) +
  labs(x = "Average CpG % methylation of gene", y = "Abs(DeSeq2 Log2FoldChange)", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

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

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->

### 0.5.4 Read counts/vsd transformed counts

RNA metadata:

``` r
meta_RNA <- read.csv("../data_RNA/LCM_RNA_metadata.csv") %>%
            dplyr::arrange(Sample) %>%
            mutate(across(c(Tissue, Fragment, Section_Date, LCM_Date), factor)) # Set variables as factors 

meta_RNA$Tissue <- factor(meta_RNA$Tissue, levels = c("OralEpi","Aboral")) #we want OralEpi to be the baseline

oral_samples <- meta_RNA$Sample[meta_RNA$Tissue == "OralEpi"]
aboral_samples <- meta_RNA$Sample[meta_RNA$Tissue == "Aboral"]
```

``` r
vsd <- read.csv("../output_RNA/differential_expression/vsd_expression_matrix.csv", row.names = 1)
vsd$vst_mean <- rowMeans(vsd)
vsd$gene_id <- rownames(vsd)

vsd <- vsd %>%select(gene_id,everything())
vsd <- vsd %>% rowwise() %>%
  mutate(Oral = mean(c_across(all_of(oral_samples)), na.rm = TRUE),
         Aboral = mean(c_across(all_of(aboral_samples)), na.rm = TRUE)) %>% 
 # mutate(Oral = mean(c_across(meta_RNA$Sample[meta_RNA$Tissue=="OralEpi"]), na.rm = TRUE)) %>%
  #mutate(Aboral = mean(c_across(meta_RNA$Sample[meta_RNA$Tissue=="Aboral"]), na.rm = TRUE)) %>%
  ungroup()
  
vsd_long <- vsd %>% pivot_longer(cols = c(Oral,Aboral),
                                      names_to = "Tissue",
                                      values_to = "tissue_vst_mean")
```

``` r
plot_data <- merge(percent_meth_long, vsd_long, by = c("gene_id", "Tissue"))

ggplot(plot_data, aes(y = tissue_vst_mean, x = tissue_percent_meth, color=Tissue)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2"))+
  labs(x = "Average CpG % methylation of gene", y = "Mean counts, all samples", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
ggplot(plot_data, aes(y = vst_mean, x = percent_meth_ALL)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "lm") + stat_poly_eq(use_label("eq", "R2"))+
  labs(x = "Average CpG % methylation of gene", y = "Mean counts, all samples", 
       title = "Gene Methylation vs Expression") +
  theme_minimal()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](09-MethylKit-ConvFilt_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
model <- lm(tissue_vst_mean ~ tissue_percent_meth * Tissue, data = plot_data)

summary(model)
```

    ## 
    ## Call:
    ## lm(formula = tissue_vst_mean ~ tissue_percent_meth * Tissue, 
    ##     data = plot_data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.7884 -1.7766 -0.2367  1.5714  8.7536 
    ## 
    ## Coefficients:
    ##                                Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    7.679166   0.048148 159.491   <2e-16 ***
    ## tissue_percent_meth            0.007577   0.002320   3.267   0.0011 ** 
    ## TissueOral                     0.672722   0.067824   9.919   <2e-16 ***
    ## tissue_percent_meth:TissueOral 0.002831   0.003356   0.844   0.3989    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.274 on 5140 degrees of freedom
    ## Multiple R-squared:  0.02764,    Adjusted R-squared:  0.02707 
    ## F-statistic:  48.7 on 3 and 5140 DF,  p-value: < 2.2e-16
