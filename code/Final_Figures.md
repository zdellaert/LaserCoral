Figures Plotted
================
Zoe Dellaert
2025-11-05

- [0.1 Setup](#01-setup)
  - [0.1.1 Load in Data](#011-load-in-data)
- [0.2 Figure 2](#02-figure-2)
  - [0.2.1 Panel B: PCA](#021-panel-b-pca)
  - [0.2.2 Panel C](#022-panel-c)
  - [0.2.3 Patchwork](#023-patchwork)
- [0.3 Figure 3](#03-figure-3)
- [0.4 Figure 4](#04-figure-4)
- [0.5 Figure 5](#05-figure-5)
- [0.6 Supplementary Tables](#06-supplementary-tables)

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.0     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.2     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(purrr)
library(ggplot2)
library(patchwork)
library(stringr)
library(scales)
```

    ## 
    ## Attaching package: 'scales'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
library(EnhancedVolcano)
```

    ## Loading required package: ggrepel

``` r
library(ComplexHeatmap)
```

    ## Loading required package: grid
    ## ========================================
    ## ComplexHeatmap version 2.26.1
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

``` r
library(RColorBrewer)
library(circlize)
```

    ## ========================================
    ## circlize version 0.4.17
    ## CRAN page: https://cran.r-project.org/package=circlize
    ## Github page: https://github.com/jokergoo/circlize
    ## Documentation: https://jokergoo.github.io/circlize_book/book/
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. circlize implements and enhances circular visualization
    ##   in R. Bioinformatics 2014.
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(circlize))
    ## ========================================

``` r
library(svglite)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:patchwork':
    ## 
    ##     align_plots
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(magick)
```

    ## Linking to ImageMagick 6.9.13.29
    ## Enabled features: cairo, fontconfig, freetype, heic, lcms, pango, raw, rsvg, webp
    ## Disabled features: fftw, ghostscript, x11

## 0.1 Setup

``` r
#set standard output directory for figures
outdir <- "../output_RNA/differential_expression"

# Specify colors
ann_colors = list(Tissue = c(OralEpi = "#FD8D3C" ,Aboral = "mediumpurple1"))
```

### 0.1.1 Load in Data

``` r
vsd <- read.csv(paste0(outdir,"/vsd_expression_matrix.csv"))
vsd <- vsd %>% column_to_rownames(var = "X")

# Read in metadata 
meta <- read.csv("../data_RNA/LCM_RNA_metadata.csv") %>%
            dplyr::arrange(Sample) %>%
            mutate(across(c(Tissue, Fragment, Section_Date, LCM_Date), factor)) %>% 
            mutate(Sample_Desc = paste0(Fragment,"_",Tissue))

# Gene annotation data

CDSearch <- read.delim("../references/Pocillopora_acuta_HIv2.genes.Conserved_Domain_Search_results.txt", quote = "") %>% dplyr::rename("query" = X.Query)
EggNog <- read.delim("../references/Pocillopora_acuta_HIv2.genes.EggNog_results.txt") %>% dplyr::rename("query" = X.query)
DESeq_SwissProt <- read.csv(file=paste0(outdir,"/DESeq_SwissProt_annotation_categorized.csv"))

DE_05_marker_broc <- read.csv(file=paste0(outdir,"/DE_05_markergene_annotation.csv"))
```

## 0.2 Figure 2

### 0.2.1 Panel B: PCA

``` r
pcaData <- read.csv(paste0(outdir, "/pcaData_allgenes.csv"), row.names = 1)
percentVar <- read.csv(paste0(outdir, "/pcaData_percentVar_allgenes.csv"))$percentVar
percentVar <- round(100 * percentVar)

PCA_small <- ggplot(pcaData, aes(PC1, PC2, color=Tissue)) +
  geom_point(size=2) +
  scale_color_manual(values = c("Aboral" = "mediumpurple1", "OralEpi" = "#FD8D3C"),
                     labels = c("Aboral" = "Aboral", "OralEpi" = "Oral"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw() + 
  theme(legend.position = c(0.8,0.25))
```

### 0.2.2 Panel C

``` r
ann_colors = list(Tissue = c(OralEpi = "#FD8D3C" ,Aboral = "mediumpurple1"))

colorRampPalette(c("mediumpurple1", "white", "#FD8D3C"))(5)
```

    ## [1] "#AB82FF" "#D4C0FF" "#FFFFFF" "#FEC69D" "#FD8D3C"

``` r
LFC_colors <- ifelse(DESeq_SwissProt$padj < 0.05 & DESeq_SwissProt$log2FoldChange > 1, '#FD8D3C',
                  ifelse(DESeq_SwissProt$padj < 0.05 & DESeq_SwissProt$log2FoldChange > 0 , '#FEC69D',
                         ifelse(DESeq_SwissProt$padj < 0.05 & DESeq_SwissProt$log2FoldChange < -1, '#AB82FF',
                                ifelse(DESeq_SwissProt$padj < 0.05 & DESeq_SwissProt$log2FoldChange < 0, '#D4C0FF',
                                       'grey'))))

LFC_colors[is.na(LFC_colors)] <- 'grey'
  names(LFC_colors)[LFC_colors == '#FD8D3C'] <- 'Up_Oral_LFC_sig'
  names(LFC_colors)[LFC_colors == '#FEC69D'] <- 'Up_Oral'
  names(LFC_colors)[LFC_colors == 'grey'] <- 'n.s.'
  names(LFC_colors)[LFC_colors == '#AB82FF'] <- 'Up_Aboral_LFC_sig'
  names(LFC_colors)[LFC_colors == '#D4C0FF'] <- 'Up_Aboral'
  
DESeq_SwissProt$Volcano_Label <-paste0(replace_na(stringr::str_wrap(DESeq_SwissProt$Heatmap_Label, width = 20), "Unannotated"),"\nLFC = ",round(DESeq_SwissProt$log2FoldChange,1))

top_DE_oral <- DESeq_SwissProt %>% filter(log2FoldChange > 1) %>% arrange(desc(log2FoldChange)) %>% head(5) %>% pull(query)
top_DE_aboral <- DESeq_SwissProt %>% filter(log2FoldChange < -1) %>% arrange(log2FoldChange) %>% head(5)  %>% pull(query)
label_top5 <- c(top_DE_oral,top_DE_aboral)

volcano <- EnhancedVolcano(DESeq_SwissProt, 
    lab = ifelse(DESeq_SwissProt$query %in% label_top5 , DESeq_SwissProt$Volcano_Label, ""),
    xlim = c(min(DESeq_SwissProt[["log2FoldChange"]], na.rm = TRUE) - 1.5, max(DESeq_SwissProt[["log2FoldChange"]], na.rm = TRUE) + 1.5),
  ylim = c(0, max(-log10(DESeq_SwissProt[["padj"]]), na.rm = TRUE) + .5),
    labSize = 3,
    axisLabSize = 12,
    captionLabSize = 4,
    x = "log2FoldChange", 
    y = "padj",
  shape=16,
    pCutoff = 0.05,
    colCustom = LFC_colors, colAlpha = .7,
    title = NULL, subtitle = NULL, caption = NULL,
    legendPosition = "none",
    legendLabSize = 8,
    pointSize = 1,
    legendIconSize = 3,
    drawConnectors = TRUE,
    lengthConnectors = unit(0.01, "npc"),
    directionConnectors = "both",
    border = 'full',
    max.overlaps = Inf,
    maxoverlapsConnectors = NULL,
    arrowheads =TRUE,
    min.segment.length = 0,
    gridlines.major = FALSE,
  gridlines.minor = FALSE,
  raster = TRUE
)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the EnhancedVolcano package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## ℹ The deprecated feature was likely used in the EnhancedVolcano package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

### 0.2.3 Patchwork

``` r
# Load images 
img_1A_crop <- ggdraw() + draw_image("images/1A_crop.png") 
img_1A_crop_cut <- ggdraw() + draw_image("images/1A_crop_cut.png")
  

top_row <- (free(img_1A_crop) + free(img_1A_crop_cut) + PCA_small) +
  plot_layout(widths = c(0.225, 0.225, 0.235)) 
final_fig <- 
  top_row / free(volcano) +
  plot_annotation(tag_levels = list(c("A1", "A2","B","C"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5),
        plot.tag.position  = c(.02, .9)) 

svglite("../output_RNA/differential_expression/composite_figs/Figure2.svg",
        width = 10, height = 7)
print(final_fig)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 0.3 Figure 3

``` r
library(ggdendro)
library(patchwork)
dist_mat <- slot(Wang_clusters_wardD2_Oral, "clusters_dist")[["BMA"]]
hc <- hclust(dist_mat, method = "ward.D2")
dend <- as.dendrogram(hc)
dend <- dendextend::rotate(dend,c(2,3,1,4,5,10,11,13,14,12,7,8,9,6))

dend_data <- ggdendro::dendro_data(dend, type = "rectangle")

dend_labels_oral <- dend_data$labels %>%
    mutate(
    cluster = str_extract(label, "\\d+"),
    Cluster.name = str_replace(label,".*?_.*?_",""),
    Cluster.term = str_replace(label,"^\\d+_",""),
    Cluster.term = str_replace(Cluster.term,"_.*$","")
  ) %>%
  mutate(Cluster.name = paste0("Cluster ", cluster, ": ", Cluster.name)) 

bar_plot_oral <- clustered_allDEGs_enrichedGO %>%  
filter(Tissue == "Oral") %>%
  arrange(as.numeric(GO.cluster)) %>%
  mutate(Cluster_label = factor(Cluster.name, levels = rev(dend_labels_oral$Cluster.name))) %>%
  ggplot(aes(x = Cluster_label)) +
  stat_count(width = .8, fill="white", color="black") +
  coord_flip() +
  theme_classic(base_size = 11) +
  scale_y_continuous(
    breaks = scales::breaks_width(2),
    expand = c(0, 0),
    limits = c(0, NA)
  ) + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, lineheight = 0.7,color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none",
    plot.subtitle = element_text(face = "bold", size = 9, hjust = 0.5),
    panel.spacing.y = unit(0.15, "lines")
  ) + labs(y = "Terms in Cluster")

dend_plot_oral <- ggplot() +
  geom_segment(
    data = dend_data$segments,
    aes(x = y, y = x, xend = yend, yend = xend),
    linewidth = 0.6,
    color = "#2C3E50"
  ) +
  scale_y_reverse(
    breaks = seq_along(dend_labels_oral$Cluster.name),
    expand = c(0.02, 0.02)
  ) +
  scale_x_reverse(expand = c(0.02, 0)) +
  theme_void()

# Combine plots with better proportions
final_oral_dend <- wrap_elements(
  dend_plot_oral + bar_plot_oral + plot_layout(widths = c(0.8, 3.5)) + theme(plot.margin = margin(0, 0, 0, 0)))
```

``` r
dist_mat <- slot(Wang_clusters_wardD2_Aboral, "clusters_dist")[["BMA"]]
hc <- hclust(dist_mat, method = "ward.D2")
dend <- as.dendrogram(hc)
dend <- dendextend::rotate(dend,c(1:8,11,12,9:10,13:16))

dend_data <- ggdendro::dendro_data(dend, type = "rectangle")

dend_labels_aboral <- dend_data$labels %>%
    mutate(
    cluster = str_extract(label, "\\d+"),
    Cluster.name = str_replace(label,".*?_.*?_",""),
    Cluster.term = str_replace(label,"^\\d+_",""),
    Cluster.term = str_replace(Cluster.term,"_.*$","")
  ) %>%
  mutate(Cluster.name = paste0("Cluster ", cluster, ": ", Cluster.name)) 

bar_plot_aboral <- clustered_allDEGs_enrichedGO %>%  
filter(Tissue == "Aboral") %>%
  arrange(as.numeric(GO.cluster)) %>%
  mutate(Cluster_label = factor(Cluster.name, levels = rev(dend_labels_aboral$Cluster.name))) %>%
  ggplot(aes(x = Cluster_label)) +
  stat_count(width = .8, fill="white", color="black") +
  coord_flip() +
  theme_classic(base_size = 11) +
  scale_y_continuous(
    breaks = scales::breaks_width(2),
    expand = c(0, 0),
    limits = c(0, NA)
  ) + 
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, lineheight = 0.7,color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "white", color = "black"),
    legend.position = "none",
    plot.subtitle = element_text(face = "bold", size = 9, hjust = 0.5),
    panel.spacing.y = unit(0.15, "lines")
  ) + labs(y = "Terms in Cluster")


dend_plot_aboral <- ggplot() +
  geom_segment(
    data = dend_data$segments,
    aes(x = y, y = x, xend = yend, yend = xend),
    linewidth = 0.6,
    color = "#2C3E50"
  ) +
  scale_y_reverse(
    breaks = seq_along(dend_labels_aboral$Cluster.name),
    expand = c(0.02, 0.02)
  ) +
  scale_x_reverse(expand = c(0.02, 0)) +
  theme_void()

# Combine plots with better proportions
final_aboral_dend <- wrap_elements(
  dend_plot_aboral + bar_plot_aboral + plot_layout(widths = c(0.8, 3.5)) + theme(plot.margin = margin(0, 0, 0, 0)))
```

``` r
top_Oral_per_Cluster <- clustered_Oral_DEGs_enrichedGO %>% group_by(GO.cluster) %>% arrange(Oral_elim.pvalue) %>% slice(1) %>% ungroup() %>% dplyr::select(GO.cluster,GO.ID,term) %>% rename("top_sig_GO" = `GO.ID`,"top_sig_term" = `term`)

plotting_oral <- left_join(clustered_Oral_DEGs_enrichedGO,top_Oral_per_Cluster) %>%
                        arrange(match(GO.cluster,dend_labels_oral$cluster)) %>%
                        mutate(is_top = GO.ID == top_sig_GO,
                               pvalue = Oral_elim.pvalue,
                                term_wrapped = str_wrap(`term`, width = 25),
                                 cluster_label = factor(str_wrap(Cluster.name,
                                        width = 20), levels = unique(str_wrap(Cluster.name,
                                        width = 20)))) %>%
  add_count(GO.cluster, name = "n_terms")

top_Aboral_per_Cluster <- clustered_Aboral_DEGs_enrichedGO %>% group_by(GO.cluster) %>% arrange(Aboral_elim.pvalue) %>% slice(1) %>% ungroup() %>% dplyr::select(GO.cluster,GO.ID,term) %>% rename("top_sig_GO" = `GO.ID`,"top_sig_term" = `term`)

plotting_aboral <- left_join(clustered_Aboral_DEGs_enrichedGO,top_Aboral_per_Cluster)%>%
                      arrange(match(GO.cluster,dend_labels_aboral$cluster)) %>%
                        mutate(is_top = GO.ID == top_sig_GO,
                               pvalue = Aboral_elim.pvalue,
                                term_wrapped = str_wrap(`term`, width = 25),
                                 cluster_label = factor(str_wrap(Cluster.name,
                                        width = 20), levels = unique(str_wrap(Cluster.name,
                                        width = 20)))) %>%
  add_count(GO.cluster, name = "n_terms")
  

combined_plotting <- rbind(plotting_oral %>%
                                        dplyr::select(-contains("Oral")),
                                      plotting_aboral %>%
                                        dplyr::select(-contains("Aboral")))
```

``` r
p_oral <-ggplot(plotting_oral, 
       aes(x = `Oral_elim.-log10_pvalue`, 
           y = reorder(GO.ID, `Oral_elim.-log10_pvalue`),
           size = `Oral_elim.-log10_pvalue`)) +
  geom_segment(aes(xend = 0, yend = reorder(GO.ID, `Oral_elim.-log10_pvalue`)),
               alpha = 0.8, linewidth = 0.6, color = "#FD8D3C") +
  geom_point(alpha = 0.8, color = "#FD8D3C") +
    facet_grid(cluster_label ~ ., 
               scales = "free_y", space = "free_y") +
    labs(x = "-log10(p-value)",
         y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 9),
          strip.text.y = element_text(size = 10, angle = 0, hjust=0,face = "bold"),
          legend.position = "none",
          panel.spacing = unit(.15, "lines"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_size_continuous(range = c(2, 5)) + 
 # scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +  # adds 10% padding on right
  coord_cartesian(clip = "off")

p_aboral <- ggplot(plotting_aboral, 
         aes(x = `Aboral_elim.-log10_pvalue`, 
             y = reorder(GO.ID, `Aboral_elim.-log10_pvalue`),
             color = factor(GO.cluster),
             size = `Aboral_elim.-log10_pvalue`)) +
    geom_segment(aes(xend = 0, yend = reorder(GO.ID, `Aboral_elim.-log10_pvalue`)),
               alpha = 0.8, linewidth = 0.6, color = "#756BB1") +
  geom_point(alpha = 0.8, color = "#756BB1") +
      facet_grid(cluster_label ~ ., 
                scales = "free_y", space = "free_y") +
  labs(x = "-log10(p-value)",
         y = NULL) +
   theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 9),
          strip.text.y = element_text(size = 10, angle = 0, hjust=0,face = "bold"),
          legend.position = "none",
          panel.spacing = unit(.15, "lines"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_size_continuous(range = c(2, 5)) + 
  coord_cartesian(clip = "off")

layout <- "
ABFF
CDFF
EEFF
EEFF
EEFF
"

final_fig <- 
  dend_plot_oral + bar_plot_oral + 
  dend_plot_aboral + bar_plot_aboral +
  free(p_oral) + free(p_aboral) + 
  plot_layout(design = layout, widths = c(1.5, 8,8,8), axes = "collect_x") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", hjust = 0, vjust = -0.5),
        plot.tag.position = c(0, 1),
        margins=margin(10.5,5.5,5.5,5.5, unit = "pt")) 

final_fig <- final_fig + plot_annotation(tag_levels = list(c("A: Oral Epidermis Tissue", "","B: Aboral Tissue", "","C: Oral Epidermis Tissue", "D: Aboral Tissue")))

ggsave("../output_RNA/differential_expression/composite_figs/Figure3.svg", final_fig, width = 12, height = 16)
```

## 0.4 Figure 4

``` r
Biomin_all_expressed <- read.csv(paste0(outdir,"/DESeq_biomin_annotation.csv")) %>% left_join(DESeq_SwissProt %>% select(query,Heatmap_Label)) %>% select(-def_short)
```

    ## Joining with `by = join_by(query)`

``` r
write.csv(Biomin_all_expressed %>% select(-X) %>% select(Heatmap_Label,query,everything()), file=paste0(outdir,"/DESeq_biomin_annotation_heatmap_labels.csv"),row.names=FALSE)

Biomin_DE <- Biomin_all_expressed %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
#choose the genes that will be plotted for the heatmap
select <- Biomin_DE$query

#gather expression info for those genes (assay(vsd) and get the z-scores)
z_scores <- t(scale(t(vsd[select, ])))
z_scores_df <- z_scores %>% as.data.frame() %>% rownames_to_column(var = "query")

Biomin_DE_Zscore <- Biomin_DE %>% left_join(z_scores_df) %>%
  rename_at(vars(meta$Sample),~ meta$Sample_Desc)
```

    ## Joining with `by = join_by(query)`

``` r
write.csv(Biomin_DE_Zscore%>% select(-X) %>% select(Heatmap_Label,query,everything()), file=paste0(outdir,"/DE05_biomin_zscore.csv"),row.names = FALSE)

# row annotations: Classification (type of biomin gene), category (genes i also called non-biomin), logfoldchange, padj
annotation_row <- Biomin_DE %>%
  select(query, Classification, log2FoldChange, padj) %>%
  tibble::column_to_rownames("query")

#row annotation: which genes are DE
annotation_row$DE <- abs(annotation_row$log2FoldChange) > 1 & annotation_row$padj < 0.05

# column annotations: tissue type
annotation_col_df <- meta %>% select(Tissue)

# colors
ann_colors = list(Tissue = c(OralEpi = "#FD8D3C" ,Aboral = "mediumpurple1"))
expr_colors <- colorRamp2(c(-2, 0, 2), c("#0571b0", "white", "#ca0020"))

Dark_Pal <- brewer.pal(length(unique(annotation_row$Classification))+2, "Dark2")

#remove purple and orange colors
Dark_Pal <- Dark_Pal[-c(2,3)]

class_colors <- setNames(
  Dark_Pal,
  unique(annotation_row$Classification)
)

# create row annotation object
row_ha <- rowAnnotation(
  Classification = annotation_row$Classification,
  col = list(
    Classification = class_colors
  ),
  annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
   annotation_legend_param = list(
    Classification = list(
      title = "Gene Classification",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  ),
  width = unit(18, "mm"),
  gap = unit(1, "mm")
)

# column annotation object
col_ha <- HeatmapAnnotation(
  Tissue = annotation_col_df$Tissue,
  col = list(Tissue = ann_colors$Tissue),
  annotation_name_gp = gpar(fontsize = 0, fontface = "bold"),
  annotation_legend_param = list(
    Tissue = list(title = "Tissue", ncol = 1),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ), height = unit(5, "mm"),
  gap = unit(1, "mm")
)

# labels
labels_all <- DESeq_SwissProt %>% 
  mutate(is_DE = abs(log2FoldChange) > 1 & padj < 0.05) %>%
  select(query, Heatmap_Label,is_DE) %>%
  mutate_all(~ ifelse(is.na(.), "", .)) #%>%
 # mutate(Heatmap_Label = ifelse(is_DE,Heatmap_Label, "")) #only label DE genes

row_labels <- labels_all$Heatmap_Label[match(rownames(z_scores), labels_all$query)]

# Perform hierarchical clustering on all samples and genes
row_hclust <- hclust(dist(z_scores))

# cut into 2 main clusters (oral vs aboral)
row_clusters <- cutree(row_hclust, k = 2)

# Combine expression-based grouping with supervised annotation
split_factor <- stringr::str_wrap(annotation_row$Classification, width = 18)
```

``` r
ht <- Heatmap(
  z_scores,
  name = "Expression\nZ-score",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_parent_dend_line = FALSE,
  row_km = 2,
  row_split = split_factor,
  column_split = annotation_col_df$Tissue,
  col = expr_colors,
  row_gap = unit(c(1,1,1,1,2,1,1,1,1), "mm"),
  #row_gap = unit(2, "mm"),
  column_gap = unit(1.5, "mm"),
  left_annotation = row_ha,
  row_title = NULL,
  top_annotation = col_ha,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6.25, fontface = "italic"),
  row_labels = row_labels,
  show_column_names = FALSE,
  column_dend_height = unit(5, "mm"),
  border = FALSE,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Expression\n(Z-Score)",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(3, "cm"),
    direction = "horizontal"
  ),
  rect_gp = gpar(col = "grey70", lwd = 0.5),
  use_raster = TRUE,
  raster_quality = 8,
  height = unit(6.25, "in")
)

svglite("../output_RNA/differential_expression/composite_figs/Figure4.svg", 
        width = 6, height = 8,
        system_fonts = list(sans = "Arial"))
draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 0.5 Figure 5

``` r
#get list of all differentially expressed marker genes
markers <- DE_05_marker_broc %>% select(query, Standardized_Name_spA) %>% distinct()

# get VST expression values
vsd_mat <- vsd %>% 
  tibble::rownames_to_column("gene")

# make VST expression value df into long format, and add marker info, filter to only keep expressed marker genes

mat_long <- vsd_mat %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "vst") %>%
  left_join(meta %>% select(Sample, Tissue), by = "Sample") %>%
  left_join(markers, by = c("gene" = "query")) %>%
  filter(!is.na(Standardized_Name_spA))
```

    ## Warning in left_join(., markers, by = c(gene = "query")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 9601 of `x` matches multiple rows in `y`.
    ## ℹ Row 43 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
# gene mean by tissue
plot_df <- mat_long %>%
  group_by(Standardized_Name_spA, gene, Tissue) %>%
  summarise(mean_vst = mean(vst, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Tissue, values_from = mean_vst, values_fill = 0) %>%
  mutate(Standardized_Name_spA = str_to_title(str_replace_all(Standardized_Name_spA,"_","\n")))
```

``` r
tissue_cols <- setdiff(colnames(plot_df), c("Standardized_Name_spA", "gene"))

# outer = OralEpi (t1), inner = Aboral (t2)
t1 <- tissue_cols[tissue_cols=="OralEpi"]
t2 <- tissue_cols[tissue_cols=="Aboral"]

# Sort sectors by tissue dominance for better visual flow
sector_totals <- plot_df %>%
  group_by(Standardized_Name_spA) %>%
  summarise(
    n_t1 = sum(.data[[t1]] > .data[[t2]], na.rm = TRUE),
    n_t2 = sum(.data[[t2]] > .data[[t1]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dom_group = if_else(n_t1 > n_t2, "t1_dominant", "t2_dominant"),
    sort_value = if_else(dom_group == "t1_dominant", -n_t1, n_t2)  # negative makes n_t1 descending
  ) %>%
  arrange(
    factor(dom_group, levels = c("t1_dominant", "t2_dominant")),
    sort_value
  )

plot_df <- plot_df %>%
  mutate(Standardized_Name_spA = factor(Standardized_Name_spA, 
                                       levels = sector_totals$Standardized_Name_spA))


gene_counts <- table(plot_df$Standardized_Name_spA)
xlims <- cbind(rep(0, length(gene_counts)), as.numeric(gene_counts))
maxval <- max(c(plot_df[[t1]], plot_df[[t2]]), na.rm = TRUE)
maxval <- signif(maxval,0)
```

``` r
#png("../output_RNA/differential_expression/DE_05_brocmarker_circ_colored.png", width = 3000, height = 3000, res = 300)
svglite("../output_RNA/differential_expression/composite_figs/Figure5.svg", width = 10, height = 10)

circos.clear()
circos.par(start.degree = 180,
           gap.after = setNames(c(rep(3,length(gene_counts)-1),8), names(gene_counts)),
           cell.padding = c(0.01, 0, 0.01, 0),
           canvas.xlim = c(-1, 1.1))
circos.initialize(factors = names(gene_counts), xlim = xlims)

# Sector labels with improved styling
circos.track(
  ylim = c(0, 1), 
  track.height = 0.08, 
  bg.border = NA,
  panel.fun = function(x, y) {
    sector <- CELL_META$sector.index
    # Calculate text size based on sector name length
    text_cex <- min(1.25, 30/nchar(sector))
    circos.text(CELL_META$xcenter, 1.85, sector,
                facing = "downward", 
                cex = text_cex,
                font = 2) 
  }
)
```

    ## Note: 1 point is out of plotting region in sector 'Neuron', track '1'.

    ## Note: 1 point is out of plotting region in sector 'Epidermis', track
    ## '1'.

    ## Note: 1 point is out of plotting region in sector 'Gland', track '1'.

    ## Note: 1 point is out of plotting region in sector 'Cnidocyte', track
    ## '1'.

    ## Note: 1 point is out of plotting region in sector 'Immune', track '1'.

    ## Note: 1 point is out of plotting region in sector 'Germline', track
    ## '1'.

    ## Note: 1 point is out of plotting region in sector 'Unidentified', track
    ## '1'.

    ## Note: 1 point is out of plotting region in sector 'Gastrodermis
    ## Symbiont Containing', track '1'.

    ## Note: 1 point is out of plotting region in sector 'Digestive
    ## Filaments', track '1'.

    ## Note: 1 point is out of plotting region in sector 'Calicoblast', track
    ## '1'.

    ## Note: 1 point is out of plotting region in sector 'Gastrodermis', track
    ## '1'.

``` r
circos.trackPlotRegion(
  ylim = c(-maxval, maxval),  # negative to positive
  track.height = 0.5,
  bg.border = "gray90",
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    genes <- plot_df %>% filter(Standardized_Name_spA == sector) %>%
                      mutate(
              dom_group = if_else(.data[[t1]] > .data[[t2]], "t1_dominant", "t2_dominant")
            ) %>%
            arrange(
              factor(dom_group, levels = c("t1_dominant", "t2_dominant")),
              desc(if_else(dom_group == "t1_dominant", .data[[t1]], -.data[[t2]]))
            )

    # Add horizontal grid lines at standardized levels
    grid_levels <- c(maxval * -0.25, maxval * -0.5, maxval * -0.75, maxval * 0.25, maxval * 0.5, maxval * 0.75)
    for (level in grid_levels) {
      circos.lines(c(0, CELL_META$xlim[2]), c(level, level),
                   col = "gray90", lwd = 0.8)
    }
    
    for (i in seq_len(nrow(genes))) {
      t1_val <- genes[[t1]][i]
      t2_val <- genes[[t2]][i]
      
      # Determine colors based on dominance
      if (t1_val > t2_val) {
        t1_col <- "#FD8D3C"  # Bright color for dominant
        t2_col <- "#E6D3F7"  # Muted color for non-dominant
      } else {
        t1_col <- "#FED8A6"  # Muted color for non-dominant
        t2_col <- "mediumpurple1"  # Bright color for dominant
      }
      
      # Tissue 1 goes up (positive)
      circos.rect(i-1, 0, i, t1_val, col = t1_col, border = "white", lwd = 0.5)
      # Tissue 2 goes down (negative) 
      circos.rect(i-1, 0, i, -t2_val, col = t2_col, border = "white", lwd = 0.5)
    }
    
        # Add y-axis labels on top-most sector
    if (CELL_META$sector.numeric.index == 1) {
      circos.yaxis(side = "left", at = c(-maxval, -maxval/2, 0, maxval/2, maxval),
                   labels = c(sprintf("%.0f", -maxval), sprintf("%.0f", -maxval/2),"0", sprintf("%.0f", maxval/2), sprintf("%.0f", maxval)),labels.niceFacing = FALSE,
                   tick.length = 0.15, labels.cex = 0.8, lwd = 1)
    }
  }
)

# Enhanced legend with better positioning and styling
par(xpd = TRUE)
legend("center", 
       legend = paste(c("Outer:", "Inner:"), c(t1, t2)), 
       fill = c("#FD8D3C", "mediumpurple1"),
       border = "black",
       bty = "o",  # Draw box around legend
       cex = 1.0,
       title = "Tissue Type",
       title.cex = 1.1)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## 0.6 Supplementary Tables

``` r
library(openxlsx)

wb <- createWorkbook()

# Table S1
table_s1 <- read_csv("../data_RNA/Table_S1_LCM_RNA_metadata_bydissection.csv")
```

    ## Rows: 10 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (3): Fragment, Tissue Type, Mean ± Standard Deviation per Dissection (µm²)
    ## dbl (1): Number of Dissections
    ## num (1): Total Tissue Area (µm²)
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
addWorksheet(wb, "Table S1")
writeData(wb, "Table S1", table_s1)

# Table S2
table_s2 <- read_csv("../output_RNA/differential_expression/DESeq_results.csv") %>% rename("query"=1)
```

    ## New names:
    ## Rows: 14464 Columns: 6
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (5): baseMean, log2FoldChange, lfcSE, pvalue, padj
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
addWorksheet(wb, "Table S2")
writeData(wb, "Table S2", table_s2)

# Table S3 (combined)
table_s3_aboral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/DE_05_Aboral_cluster_heatmap_Wang_wardD2.tsv") %>%
  rename_with(~gsub("^Aboral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Aboral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 123 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Aboral_elim.genes_frequency, Aboral_elim.S...
    ## dbl (4): GO.cluster, IC, Aboral_elim.pvalue, Aboral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s3_oral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/DE_05_Oral_cluster_heatmap_Wang_wardD2.tsv") %>% 
  rename_with(~gsub("^Oral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Oral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 63 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Oral_elim.genes_frequency, Oral_elim.Signi...
    ## dbl (4): GO.cluster, IC, Oral_elim.pvalue, Oral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s3 <- bind_rows(table_s3_aboral, table_s3_oral)
addWorksheet(wb, "Table S3")
writeData(wb, "Table S3", table_s3)

# Table S4 (combined)
table_s4_aboral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/MF_DE_05_Aboral_cluster_heatmap_Wang_wardD2.tsv") %>% 
  rename_with(~gsub("^Aboral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Aboral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 27 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Aboral_elim.genes_frequency, Aboral_elim.S...
    ## dbl (4): GO.cluster, IC, Aboral_elim.pvalue, Aboral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s4_oral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/MF_DE_05_Oral_cluster_heatmap_Wang_wardD2.tsv") %>% 
  rename_with(~gsub("^Oral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Oral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 28 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Oral_elim.genes_frequency, Oral_elim.Signi...
    ## dbl (4): GO.cluster, IC, Oral_elim.pvalue, Oral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s4 <- bind_rows(table_s4_aboral, table_s4_oral)
addWorksheet(wb, "Table S4")
writeData(wb, "Table S4", table_s4)

# Table S5 (combined)
table_s5_aboral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/CC_DE_05_Aboral_cluster_heatmap_Wang_wardD2.tsv") %>% 
  rename_with(~gsub("^Aboral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Aboral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 26 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Aboral_elim.genes_frequency, Aboral_elim.S...
    ## dbl (4): GO.cluster, IC, Aboral_elim.pvalue, Aboral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s5_oral <- read_tsv("../output_RNA/differential_expression/semantic-enrichment/CC_DE_05_Oral_cluster_heatmap_Wang_wardD2.tsv") %>% 
  rename_with(~gsub("^Oral_", "", .)) %>%
  mutate(Tissue_Upregulation = "Oral") %>%
  select(Tissue_Upregulation,everything())
```

    ## Rows: 19 Columns: 10
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): GO.ID, term, definition, Oral_elim.genes_frequency, Oral_elim.Signi...
    ## dbl (4): GO.cluster, IC, Oral_elim.pvalue, Oral_elim.-log10_pvalue
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
table_s5 <- bind_rows(table_s5_aboral, table_s5_oral)
addWorksheet(wb, "Table S5")
writeData(wb, "Table S5", table_s5)

## Table S6
table_s6 <- read_csv("../output_RNA/differential_expression/DESeq_biomin_annotation_heatmap_labels.csv")
```

    ## Rows: 182 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): Heatmap_Label, query, List, definition, Classification, Reference, ...
    ## dbl (5): baseMean, log2FoldChange, lfcSE, pvalue, padj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
addWorksheet(wb, "Table S6")
writeData(wb, "Table S6", table_s6)

## Table S7
table_s7 <- read_csv("../output_RNA/differential_expression/DE05_biomin_zscore.csv")
```

    ## Rows: 72 Columns: 22
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (7): Heatmap_Label, query, List, definition, Classification, Reference,...
    ## dbl (15): baseMean, log2FoldChange, lfcSE, pvalue, padj, C_OralEpi, C_Aboral...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
addWorksheet(wb, "Table S7")
writeData(wb, "Table S7", table_s7)

## Table S8
table_s8 <- read_csv("../output_RNA/differential_expression/DESeq_markergene_annotation.csv") %>% select(-1) 
```

    ## New names:
    ## Rows: 238 Columns: 9
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (3): query, Spis_CellType_Full, Spis_CellType dbl (6): ...1, baseMean,
    ## log2FoldChange, lfcSE, pvalue, padj
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
addWorksheet(wb, "Table S8")
writeData(wb, "Table S8", table_s8)

saveWorkbook(wb, "../Supplementary_Tables.xlsx", overwrite = TRUE)
```
