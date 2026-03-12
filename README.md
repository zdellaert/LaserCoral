# LaserCoral
Data and methods for Chapter Three of my dissertation, Laser Capture Microdissection RNAseq analysis

## Repository organization

**Sample metadata**
- [Metadata for all samples](https://github.com/zdellaert/LaserCoral/blob/main/LCM_RNA_metadata.csv)

**RNA-seq Data QC and Processing**
- [Bioinformatic pipeline](https://github.com/zdellaert/LaserCoral/blob/main/code/RNA-seq-bioinf.md)   
- QC of [Raw](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/raw_qc) and [Trimmed](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/trimmed_oligo_qc) reads
- [Raw data: Sequence upload to SRA, BioProject PRJNA1209584](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1209584)
- [Genome mapping rate per sample](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/hisat2/mapped_reads_counts_Pacuta)
- [Gene count matrix](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/stringtie-GeneExt/LCM_RNA_gene_count_matrix.csv)

**DeSeq2 differential expression**
- [DeSeq2 Differential Expression Analysis script](https://github.com/zdellaert/LaserCoral/blob/main/code/02-DE.md)
- [GO Enrichment Analysis](https://github.com/zdellaert/LaserCoral/blob/main/code/06-Semantic-Enrichment.Rmd)
- [Biomineralization + Marker Gene Analysis (and all final figures)](https://github.com/zdellaert/LaserCoral/blob/main/code/Final_Figures.md)
- [Outputs](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/differential_expression)
  - [GO Enrichment Outputs](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/differential_expression/semantic-enrichment)

**Biomineralization-related genes ([Scucchia et al., 2021](https://doi.org/10.1111/gcb.15812); [Scucchia et al., 2021](https://doi.org/10.1098/rspb.2021.0328))**
- [Excel file](https://github.com/zdellaert/LaserCoral/tree/main/references/Biomineralization_Toolkit_FScucchia/Biomineralization_Toolkit_FScucchia.xlsx)
- [Fasta file](https://github.com/zdellaert/LaserCoral/tree/main/references/Biomineralization_Toolkit_FScucchia/Biomineralization_Toolkit_FScucchia.fasta)
- [*P. acuta* orthologs](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/marker_genes/Pacuta_Biomin_Spis_ortholog.csv)

**Stylophora pistillata scRNA-seq Marker Genes ([Levy et al., 2021](10.1016/j.cell.2021.04.005)**
- [Marker Gene File](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/marker_genes/Spis_mc_coral_markers_list_cell_type.tsv) downloaded from [here](https://github.com/sebepedroslab/Stylophora_single_cell_atlas/blob/master/clustering_coral/scdb/Spis_mc_coral_markers_list_cell_type.tsv)
- [*P. acuta* orthologs](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/marker_genes/Pacuta_Spis_Markers_pairs.csv)

***Pocillopora acuta* species identification:**
- [FASTA and genious files of *Pocillopora* species tress](https://github.com/zdellaert/LaserCoral/tree/main/data_speciesID)

**Reference Genome and annotations:**   
- Version 2 from: http://cyanophora.rutgers.edu/Pocillopora_acuta/

**NCBI Sequence uploads:**
- Bioproject: PRJNA1209584  