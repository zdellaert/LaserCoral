Enrichment
================
Zoe Dellaert
2024-12-09

- [0.1 Gene Ontology Analysis analysis of LCM RNA
  Data](#01-gene-ontology-analysis-analysis-of-lcm-rna-data)
- [0.2 Managing Packages Using Renv](#02-managing-packages-using-renv)
- [0.3 Load packages](#03-load-packages)
- [0.4 Description of pipeline](#04-description-of-pipeline)
- [0.5 Load in reference files and differential expression
  data](#05-load-in-reference-files-and-differential-expression-data)
- [0.6 Create gene_to_go.tab file for running
  GO_MWU](#06-create-gene_to_gotab-file-for-running-go_mwu)
- [0.7 Create the GO_MWU input csv files for the binary
  analysis](#07-create-the-go_mwu-input-csv-files-for-the-binary-analysis)
- [0.8 Create the GO_MWU input csv files for the continuous LFC
  analysis](#08-create-the-go_mwu-input-csv-files-for-the-continuous-lfc-analysis)
- [0.9 Run GO_MWU](#09-run-go_mwu)
  - [0.9.1 Up Aboral Fisher Exact
    Test](#091-up-aboral-fisher-exact-test)
  - [0.9.2 Oral Epidermis tissue fisher exact
    test](#092-oral-epidermis-tissue-fisher-exact-test)
  - [0.9.3 Differential Expression using
    LFC](#093-differential-expression-using-lfc)
- [0.10 Custom plots](#010-custom-plots)

## 0.1 Gene Ontology Analysis analysis of LCM RNA Data

## 0.2 Managing Packages Using Renv

To run this code in my project using the renv environment, run the
following lines of code

``` r
install.packages("renv") #install the package on the new computer (may not be necessary if renv bootstraps itself as expected)
renv::restore() #reinstall all the package versions in the renv lockfile
```

## 0.3 Load packages

``` r
require("ape")
```

    ## Loading required package: ape

    ## Warning: package 'ape' was built under R version 4.3.3

``` r
require("scales")
```

    ## Loading required package: scales

``` r
require("tidyplots")
```

    ## Loading required package: tidyplots

    ## Warning: package 'tidyplots' was built under R version 4.3.3

``` r
require("tidyverse")
```

    ## Loading required package: tidyverse

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ readr::col_factor() masks scales::col_factor()
    ## ✖ purrr::discard()    masks scales::discard()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ dplyr::where()      masks ape::where()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
sessionInfo() #provides list of loaded packages and version of R.
```

    ## R version 4.3.2 (2023-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.0
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices datasets  utils     methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
    ##  [5] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
    ##  [9] ggplot2_3.5.1   tidyverse_2.0.0 tidyplots_0.2.0 scales_1.3.0   
    ## [13] ape_5.8-1      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.5        compiler_4.3.2      BiocManager_1.30.25
    ##  [4] renv_1.0.11         tidyselect_1.2.1    Rcpp_1.0.13-1      
    ##  [7] parallel_4.3.2      yaml_2.3.10         fastmap_1.2.0      
    ## [10] lattice_0.22-6      R6_2.5.1            generics_0.1.3     
    ## [13] knitr_1.48          munsell_0.5.1       tzdb_0.4.0         
    ## [16] pillar_1.9.0        rlang_1.1.4         utf8_1.2.4         
    ## [19] stringi_1.8.4       xfun_0.48           timechange_0.3.0   
    ## [22] cli_3.6.3           withr_3.0.1         magrittr_2.0.3     
    ## [25] digest_0.6.37       grid_4.3.2          rstudioapi_0.17.0  
    ## [28] hms_1.1.3           lifecycle_1.0.4     nlme_3.1-166       
    ## [31] vctrs_0.6.5         evaluate_1.0.1      glue_1.8.0         
    ## [34] fansi_1.0.6         colorspace_2.1-1    rmarkdown_2.28     
    ## [37] tools_4.3.2         pkgconfig_2.0.3     htmltools_0.5.8.1

## 0.4 Description of pipeline

I am going to perform functional enrichment of GO terms using
[GO_MWU](https://github.com/z0on/GO_MWU).

As described in the [README](https://github.com/z0on/GO_MWU) of the
package, the steps to running this analysis are as follows:

1.  Make a directory for GO_MWU files. I am creating one in this code/
    directory: LaserCoral/code/go_mwu
2.  Download the following scripts into that directory: GO_MWU.R,
    gomwu_a.pl, gomwu_b.pl, gomwu.functions.R (downloaded from the
    github, linked above)
3.  Download the GO hierarchy file, go.obo, from this website:
    <http://www.geneontology.org/GO.downloads.ontology.shtml>
4.  table of GO annotations for your sequences:

- “two-column (gene id - GO terms), tab-delimited, one line per gene,
  multiple GO terms separated by semicolon. If you have multiple lines
  per gene, use nrify_GOtable.pl to merge them. Do NOT include genes
  without GO annotations.”

5.  table of measure of interest for your sequences:

- “two columns of comma-separated values: gene id, continuous measure of
  change such as log(fold-change). To perform standard GO enrichment
  analysis based on Fisher’s exact test, use binary measure (1 or 0,
  i.e., either significant or not).”

## 0.5 Load in reference files and differential expression data

In the next chunk I am loadinf in my DESeq data. These results are
ordered by adjusted p-value. As a reminder, negative LFC = higher in
Aboral tissue, and positive LFC = higher in Oral tissue.

``` r
#load in DESeq results
DESeq <- read.csv("../output_RNA/differential_expression/DESeq_results.csv", header = TRUE) %>% dplyr::rename("query" ="X")

#make dataframes of just differentially expressed genes for each LFC direction
DE_05_Aboral <- DESeq %>% filter(padj < 0.05 & log2FoldChange > 0)
DE_05_OralEpi <- DESeq %>% filter(padj < 0.05& log2FoldChange < 0)

#load in annotation data 
EggNog <- read.delim("../references/Pocillopora_acuta_HIv2.genes.EggNog_results.txt") %>% dplyr::rename("query" = X.query) 

#filter annotation data for just expressed genes 
EggNog <- EggNog %>% filter(query %in% DESeq$query)
EggNog$GOs <- gsub("-", NA, EggNog$GOs)
EggNog$GOs <- gsub(",", ";", EggNog$GOs)

nrow(EggNog)
```

    ## [1] 9482

``` r
nrow(EggNog)/nrow(DESeq)
```

    ## [1] 0.6555586

Only 9482/14464 genes in our dataset have annotation information in this
file. That is 66%.

``` r
sum(EggNog$query %in% DE_05_Aboral$query)
```

    ## [1] 383

``` r
sum(EggNog$query %in% DE_05_Aboral$query)/nrow(DE_05_Aboral)
```

    ## [1] 0.4763682

``` r
sum(EggNog$query %in% DE_05_OralEpi$query)
```

    ## [1] 1689

``` r
sum(EggNog$query %in% DE_05_OralEpi$query)/nrow(DE_05_OralEpi)
```

    ## [1] 0.6027837

Only 383/804 genes that are significantly upregulated in the Aboral
tissue have annotation information. That is only 47% of the genes.

Only 1689/2802 genes that are significantly upregulated in the Oral
Epidermis tissue have annotation information. That is 60% of the genes.

## 0.6 Create gene_to_go.tab file for running GO_MWU

``` r
GO.terms <- EggNog %>% select(query,GOs) %>% dplyr::rename("GO.terms" = GOs)

GO.terms_only <- GO.terms %>% na.omit()
nrow(GO.terms_only)
```

    ## [1] 6580

``` r
nrow(GO.terms_only)/nrow(DESeq)
```

    ## [1] 0.4549226

``` r
write.table(GO.terms_only, "go_mwu/gene_to_go.tab",row.names = FALSE, sep = "\t", quote = FALSE)
```

Only 6580/14464 genes in our dataset have GO term information. That is
less than half: 45%.

``` r
sum(GO.terms_only$query %in% DE_05_Aboral$query)
```

    ## [1] 214

``` r
sum(GO.terms_only$query %in% DE_05_Aboral$query)/nrow(DE_05_Aboral)
```

    ## [1] 0.2661692

``` r
sum(GO.terms_only$query %in% DE_05_OralEpi$query)
```

    ## [1] 1045

``` r
sum(GO.terms_only$query %in% DE_05_OralEpi$query)/nrow(DE_05_OralEpi)
```

    ## [1] 0.3729479

Only 214/804 genes that are significantly upregulated in the Aboral
tissue have GO term information. That is only 27% of the genes.

Only 1045/2802 genes that are significantly upregulated in the Oral
Epidermis tissue have GO term information. That is 37% of the genes.

## 0.7 Create the GO_MWU input csv files for the binary analysis

``` r
### Generate vector with 1 for the 214 significant DEGs that have positive LFC and GO annotation, and 0 for the other genes that have GO annotation

UpAboral <- GO.terms_only %>% inner_join(DESeq) %>%
    mutate(DE = ifelse(padj < 0.05 & log2FoldChange > 0, 1, 0)) %>%
    dplyr::select(query, DE) 
```

    ## Joining with `by = join_by(query)`

``` r
print(sum(UpAboral$DE)) #only 214 of these have GO annotations
```

    ## [1] 214

``` r
write.csv(UpAboral, "go_mwu/UpAboral.csv", row.names = FALSE, quote = FALSE)

### Generate vector with 1 for the 1045 significant DEGs that have positive LFC and GO annotation, and 0 for the other genes that have GO annotation

UpOralEpi <- GO.terms_only %>% inner_join(DESeq) %>%
    mutate(DE = ifelse(padj < 0.05 & log2FoldChange < 0, 1, 0)) %>%
    dplyr::select(query, DE) 
```

    ## Joining with `by = join_by(query)`

``` r
print(sum(UpOralEpi$DE)) #only 1045 of these have GO annotations
```

    ## [1] 1045

``` r
write.csv(UpOralEpi, "go_mwu/UpOralEpi.csv", row.names = FALSE, quote = FALSE)
```

## 0.8 Create the GO_MWU input csv files for the continuous LFC analysis

``` r
### Generate vector with all of the genes that are expressed and have GO annotation, with their Log2FoldChanges
DE_LFC <- GO.terms_only %>% inner_join(DESeq) %>%
    dplyr::select(query, log2FoldChange) 
```

    ## Joining with `by = join_by(query)`

``` r
write.csv(DE_LFC, "go_mwu/DE_LFC.csv", row.names = FALSE, quote = FALSE)
```

## 0.9 Run GO_MWU

``` r
setwd("go_mwu")

goAnnotations="gene_to_go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")
```

### 0.9.1 Up Aboral Fisher Exact Test

``` r
setwd("go_mwu")
input="UpAboral.csv"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
#   ,Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)
```

    go.obo gene_to_go.tab UpAboral.csv BP largest=0.1 smallest=5 cutHeight=0.25

    Run parameters:

    largest GO category as fraction of all genes (largest)  : 0.1
             smallest GO category as # of genes (smallest)  : 5
                    clustering threshold (clusterCutHeight) : 0.25

    -----------------
    retrieving GO hierarchy, reformatting data...

    -------------
    go_reformat:
    Genes with GO annotations, but not listed in measure table: 1

    Terms without defined level (old ontology?..): 969
    -------------
    -------------
    go_nrify:
    13950 categories, 6284 genes; size range 5-628.4
        214 too broad
        6258 too small
        7478 remaining

    removing redundancy:

    calculating GO term similarities based on shared genes...
    5277 non-redundant GO categories of good size
    -------------

    Secondary clustering:
    calculating similarities....
    Binary classification detected; will perform Fisher's test
    331 GO terms at 10% FDR

#### 0.9.1.1 Plotting

``` r
setwd("go_mwu")
input="UpAboral.csv"

png(filename = paste0(input, "_plot.png"), width = 1200, height = 4800)
results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=0.001,  #Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment
    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
    level2=0.01, # FDR cutoff to print in regular (not italic) font.
    level3=0.001, # FDR cutoff to print in large bold font.
    txtsize=2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
    treeHeight=0.5, # height of the hierarchical clustering tree
 )
```

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "", :
    ## the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "", :
    ## the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## GO terms dispayed: 243

    ## "Good genes" accounted for:  152 out of 202 ( 75% )

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# text representation of results, with actual adjusted p-values
head(results[[1]])
```

    ##                                                          pval direction  color
    ## 8/54 proteoglycan metabolic process              1.164399e-02         0 grey50
    ## 13/126 carbohydrate derivative catabolic process 8.425550e-03         0  black
    ## 10/43 aminoglycan catabolic process              1.544289e-04         0  black
    ## 3/5 chitin catabolic process                     1.164399e-02         0 grey50
    ## 5/12 chitin metabolic process                    1.798878e-03         0  black
    ## 16/80 aminoglycan metabolic process              6.157244e-06         0  black

#### 0.9.1.2 extracting representative GOs

``` r
# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
 plot(results[[2]],cex=0.3)
 abline(h=hcut,col="red")
```

![](05-Enrichment_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
    rn=names(ct)[ct==ci]
    obs=grep("obsolete",rn)
    if(length(obs)>0) { rn=rn[-obs] }
    if (length(rn)==0) {next}
    rr=results[[1]][rn,]
    bestrr=rr[which(rr$pval==min(rr$pval)),]
    best=1
    if(nrow(bestrr)>1) {
        nns=sub(" .+","",row.names(bestrr))
        fr=c()
        for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
        best=which(fr==max(fr))
    }
    if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}
```

    ## 1

    ## 2

    ## 3

``` r
setwd("go_mwu")
mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs
```

    ##      delta.rank         pval level nseqs                  term
    ## 452           0 3.685232e-09     3    80 GO:0006022;GO:0030203
    ## 3095          0 2.353664e-06     2   140            GO:0045229
    ## 5095          0 2.165615e-11     2    13 GO:1905483;GO:1905485
    ##                                               name        p.adj
    ## 452                  aminoglycan metabolic process 6.157244e-06
    ## 3095 external encapsulating structure organization 3.184085e-04
    ## 5095          regulation of motor neuron migration 1.142579e-07

------------------------------------------------------------------------

### 0.9.2 Oral Epidermis tissue fisher exact test

``` r
setwd("go_mwu")
input="UpOralEpi.csv"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
#   ,Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)
```

    go.obo gene_to_go.tab UpOralEpi.csv BP largest=0.1 smallest=5 cutHeight=0.25

    Run parameters:

    largest GO category as fraction of all genes (largest)  : 0.1
             smallest GO category as # of genes (smallest)  : 5
                    clustering threshold (clusterCutHeight) : 0.25

    -----------------
    retrieving GO hierarchy, reformatting data...

    -------------
    go_reformat:
    Genes with GO annotations, but not listed in measure table: 1

    Terms without defined level (old ontology?..): 969
    -------------
    -------------
    go_nrify:
    13950 categories, 6284 genes; size range 5-628.4
        214 too broad
        6258 too small
        7478 remaining

    removing redundancy:

    calculating GO term similarities based on shared genes...
    5277 non-redundant GO categories of good size
    -------------

    Secondary clustering:
    calculating similarities....
    Binary classification detected; will perform Fisher's test
    0 GO terms at 10% FDR

**0 GO terms at 10% FDR**

#### 0.9.2.1 No further analysis for this LFC direction.

------------------------------------------------------------------------

### 0.9.3 Differential Expression using LFC

``` r
setwd("go_mwu")
input="DE_LFC.csv"

gomwuStats(input, goDatabase, goAnnotations, goDivision,
    perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
    largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
    smallest=5,   # a GO category should contain at least this many genes to be considered
    clusterCutHeight=0.25 # threshold for merging similar (gene-sharing) terms. See README for details.
#   ,Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
)
```

    go.obo gene_to_go.tab DE_LFC.csv BP largest=0.1 smallest=5 cutHeight=0.25

    Run parameters:

    largest GO category as fraction of all genes (largest)  : 0.1
             smallest GO category as # of genes (smallest)  : 5
                    clustering threshold (clusterCutHeight) : 0.25

    -----------------
    retrieving GO hierarchy, reformatting data...

    -------------
    go_reformat:
    Genes with GO annotations, but not listed in measure table: 1

    Terms without defined level (old ontology?..): 969
    -------------
    -------------
    go_nrify:
    13950 categories, 6284 genes; size range 5-628.4
        214 too broad
        6258 too small
        7478 remaining

    removing redundancy:

    calculating GO term similarities based on shared genes...
    5277 non-redundant GO categories of good size
    -------------

    Secondary clustering:
    calculating similarities....
    Continuous measure of interest: will perform MWU test
    602 GO terms at 10% FDR

#### 0.9.3.1 Plotting

``` r
setwd("go_mwu")
input="DE_LFC.csv"

png(filename = paste0(input, "_plot.png"), width = 1200, height = 4800)
results=gomwuPlot(input,goAnnotations,goDivision,
    absValue=1, # un-remark this if you are using log2-fold changes
    level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
    level2=0.01, # FDR cutoff to print in regular (not italic) font.
    level3=0.001, # FDR cutoff to print in large bold font.
    txtsize=1.9,    # decrease to fit more on one page
    treeHeight=0.5, # height of the hierarchical clustering tree
    colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") 
 )
```

    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "", :
    ## the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'
    ## Warning in plot.formula(c(1:top) ~ c(1:top), type = "n", axes = F, xlab = "", :
    ## the formula 'c(1:top) ~ c(1:top)' is treated as 'c(1:top) ~ 1'

    ## GO terms dispayed: 436

    ## "Good genes" accounted for:  479 out of 632 ( 76% )

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# text representation of results, with actual adjusted p-values
head(results[[1]])
```

    ##                                                                                                              pval
    ## 0/9 endonucleolytic cleavage in 5'-ETS of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA) 0.04144612
    ## 2/33 spliceosomal snRNP assembly                                                                       0.02858913
    ## 4/43 regulation of DNA recombination                                                                   0.03375225
    ## 3/20 meiotic chromosome separation                                                                     0.04562801
    ## 5/80 nucleotide-excision repair                                                                        0.02769915
    ## 3/6 DNA unwinding involved in DNA replication                                                          0.02869358
    ##                                                                                                        direction
    ## 0/9 endonucleolytic cleavage in 5'-ETS of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)         1
    ## 2/33 spliceosomal snRNP assembly                                                                               0
    ## 4/43 regulation of DNA recombination                                                                           0
    ## 3/20 meiotic chromosome separation                                                                             0
    ## 5/80 nucleotide-excision repair                                                                                0
    ## 3/6 DNA unwinding involved in DNA replication                                                                  0
    ##                                                                                                             color
    ## 0/9 endonucleolytic cleavage in 5'-ETS of tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA) lightcoral
    ## 2/33 spliceosomal snRNP assembly                                                                         skyblue2
    ## 4/43 regulation of DNA recombination                                                                     skyblue2
    ## 3/20 meiotic chromosome separation                                                                       skyblue2
    ## 5/80 nucleotide-excision repair                                                                          skyblue2
    ## 3/6 DNA unwinding involved in DNA replication                                                            skyblue2

#### 0.9.3.2 extracting representative GOs

``` r
# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
 plot(results[[2]],cex=0.1)
 abline(h=hcut,col="red")
```

![](05-Enrichment_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
    rn=names(ct)[ct==ci]
    obs=grep("obsolete",rn)
    if(length(obs)>0) { rn=rn[-obs] }
    if (length(rn)==0) {next}
    rr=results[[1]][rn,]
    bestrr=rr[which(rr$pval==min(rr$pval)),]
    best=1
    if(nrow(bestrr)>1) {
        nns=sub(" .+","",row.names(bestrr))
        fr=c()
        for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
        best=which(fr==max(fr))
    }
    if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}
```

    ## 1

    ## 2

    ## 3

    ## 4

    ## 5

    ## 6

    ## 7

    ## 8

    ## 9

    ## 10

    ## 11

    ## 12

    ## 13

    ## 14

    ## 15

``` r
setwd("go_mwu")
mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs
```

    ##      delta.rank         pval level nseqs
    ## 256        1230 8.595238e-11     2    91
    ## 452        1149 1.276802e-08     2    80
    ## 705        -996 1.599151e-05     2    61
    ## 4553       1129 2.059403e-04     2    35
    ## 4988       1168 1.774295e-07     2    65
    ## 5086       1169 6.916407e-11     2   102
    ##                                                        term
    ## 256                                              GO:0002181
    ## 452                                   GO:0006022;GO:0030203
    ## 705  GO:0006754;GO:0009206;GO:0009145;GO:0009201;GO:0009142
    ## 4553                                             GO:0097006
    ## 4988                                             GO:1903725
    ## 5086                                             GO:1905330
    ##                                                  name        p.adj
    ## 256                           cytoplasmic translation 2.267424e-07
    ## 452                     aminoglycan metabolic process 8.420510e-06
    ## 705      nucleoside triphosphate biosynthetic process 1.298019e-03
    ## 4553 regulation of plasma lipoprotein particle levels 6.280583e-03
    ## 4988     regulation of phospholipid metabolic process 4.586727e-05
    ## 5086     regulation of morphogenesis of an epithelium 2.267424e-07

## 0.10 Custom plots

``` r
MWU_BP_DE_LFC <- read.table("go_mwu/MWU_BP_DE_LFC.csv", header = TRUE)
BP_DE_LFC <- read.table("go_mwu/BP_DE_LFC.csv", header = TRUE)
```

``` r
MWU_BP_DE_LFC %>% dplyr::filter(p.adj<0.05 & delta.rank > 0) %>% 
  mutate(log10padj = -log10(p.adj)) %>%
  tidyplot(y = name, x = nseqs, fill = log10padj) %>%
  adjust_size(height = NA, width=NA) %>%
  theme_minimal_x() %>%
  add_sum_bar() %>% #adjust_colors(new_colors = colors_continuous_mako) %>%
  sort_y_axis_labels(delta.rank) %>%
  adjust_padding(right = 0.1) %>%
  add_data_labels(delta.rank, label_position = c("right")) %>%
  adjust_legend_title("$-log_[10]~(Adjusted~p-value)$") %>%
  add_title("Enriched GO terms by GO_MWU, Up in Aboral Tissue, padj < 0.05") %>%
  save_plot("../output_RNA/differential_expression/enrichment/GO_MWU_UpAboral_05.png",
            width = 8, height = 30, units="in",bg = "transparent") 
```

    ## ✔ save_plot: saved to '../output_RNA/differential_expression/enrichment/GO_MWU_UpAboral_05.png'

![](05-Enrichment_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
MWU_BP_DE_LFC %>% dplyr::filter(p.adj<0.001 & delta.rank > 0) %>% 
  mutate(log10padj = -log10(p.adj)) %>%
  tidyplot(y = name, x = nseqs, fill = log10padj) %>%
  adjust_size(height = NA, width=NA) %>%
  theme_minimal_x() %>%
  add_sum_bar() %>% #adjust_colors(new_colors = colors_continuous_mako) %>%
  sort_y_axis_labels(delta.rank) %>%
  adjust_padding(right = 0.1) %>%
  add_data_labels(delta.rank, label_position = c("right")) %>%
  adjust_legend_title("$-log_[10]~(Adjusted~p-value)$") %>%
  add_title("Enriched GO terms by GO_MWU, Up in Aboral Tissue, padj < 0.001") %>%
  save_plot("../output_RNA/differential_expression/enrichment/GO_MWU_UpAboral.png",
            width = 8, height = 6, units="in",bg = "transparent") 
```

    ## ✔ save_plot: saved to '../output_RNA/differential_expression/enrichment/GO_MWU_UpAboral.png'

![](05-Enrichment_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
MWU_BP_DE_LFC %>% dplyr::filter(p.adj<0.05 & delta.rank < 0) %>% 
  mutate(log10padj = -log10(p.adj)) %>%
  tidyplot(y = name, x = nseqs, fill = log10padj) %>%
  adjust_size(height = NA, width=NA) %>%
  theme_minimal_x() %>%
  add_sum_bar() %>% #adjust_colors(new_colors = colors_continuous_mako) %>%
  sort_y_axis_labels(-delta.rank) %>%
  adjust_padding(right = 0.1) %>%
  add_data_labels(delta.rank, label_position = c("right")) %>%
  adjust_legend_title("$-log_[10]~(Adjusted~p-value)$") %>%
  add_title("Enriched GO terms by GO_MWU, Up in Oral Epidermis Tissue, padj < 0.05")  %>%
  save_plot("../output_RNA/differential_expression/enrichment/GO_MWU_UpOralEpi.png",
            width = 8, height = 6, units="in",bg = "transparent")
```

    ## ✔ save_plot: saved to '../output_RNA/differential_expression/enrichment/GO_MWU_UpOralEpi.png'

![](05-Enrichment_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->
