# LCM RNA-seq Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 10/2/2024

- Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
- RNA extractions: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Extractions/
- RNA library prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/

## Make directory structure on Andromeda

```
mkdir /data/putnamlab/zdellaert/LaserCoral
cd /data/putnamlab/zdellaert/LaserCoral #Enter working directory

mkdir scripts #make folder for scripts
mkdir data_RNA  #make folder for raw data
mkdir output #make folder for output
```

## Check integrity of data transfer

raw data is located in `/data/putnamlab/KITT/hputnam/20241002_LaserCoral`

```
cd data_RNA
cp /data/putnamlab/KITT/hputnam/20241002_LaserCoral/20241002_URI.md5 .
cat /data/putnamlab/KITT/hputnam/20241002_LaserCoral/*gz.md5 > genohub.md5

#use diff command to see if there is a differnece between the checksums
# using -w to ignore spaces (the URI file is formatted slightly differently) 

diff -w genohub.md5 20241002_URI.md5 

# no output, so no difference between the md5s. verified by inspecting manually.
```

Data appears to have been transferred successfully from genohub.

## Symlink raw data files into data_RNA

```
ln -s /data/putnamlab/KITT/hputnam/20241002_LaserCoral/*.fastq.gz /data/putnamlab/zdellaert/LaserCoral/data_RNA
```

## QC raw files

```
cd /data/putnamlab/zdellaert/LaserCoral #Enter working directory
nano scripts/raw_qc.sh #write script for QC, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error=outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/data_RNA

# load modules needed
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make raw_qc output folder
mkdir raw_qc/

# Make an array of fastq files to trim
array1=($(ls *.fastq.gz)) 

#run fastqc on raw data
for i in ${array1[@]}; do
    fastqc ${i} -o raw_qc/
    echo "fastqc ${i} done"
done

#fastqc *.fastq.gz -o raw_qc/

#Compile MultiQC report from FastQC files
multiqc raw_qc/  #Compile MultiQC report from FastQC files 

mv multiqc_report.html raw_qc/raw_qc_multiqc_report.html
mv multiqc_data raw_qc/raw_multiqc_data

echo "Initial QC of Seq data complete." $(date)
```

```
sbatch scripts/raw_qc.sh
```

### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/raw_qc)

Okay! So there is definitely data! Yay!

And a lot of it!

And it isn't horrible quality! Yay!

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot1.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot2.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot3.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot4.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot5.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot6.png?raw=true">

But there is also a lot of duplication, overrepresentation of seqeunces, and adapter content. This isn't unexpected given the fact that this was low-input of degraded RNA extracted from fixed, LCM-d tissue. Let's see what we can do to interpret and filter from here.

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot_fastqc_adapter.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/raw_qc/multiqc_screenshots/screenshot_fastqc_overrep.png?raw=true">

**This adapter content graph is indicative that many of the inserts were much shorter than 150bp**

In the [library preparation](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/), several [oligos and adapters](https://github.com/zdellaert/ZD_Putnam_Lab_Notebook/blob/master/protocols/manualE6420_NEBNext_Low_Input_RNA_Library_Prep.pdf) are used:

1. Oligo Sequences
   1. NEBNext Template Switching Oligo: 5 ́-GCT AAT CAT TGC AAG CAG TGG TAT CAA CGC AGA GTA CAT rGrGrG-3 ́
   2. NEBNext Single Cell RT Primer Mix: 5 ́-AAG CAG TGG TAT CAA CGC AGA GTA CTT TTT TTT TTT TTT TTT TTT TTT TTT TTT TV-3 ́
   3. NEBNext Single Cell cDNA PCR Primer: 5 ́-AAG CAG TGG TAT CAA CGC AGA GT-3 ́

2. Adaptor Trimming Sequences: The NEBNext libraries for Illumina resemble TruSeq libraries and can be trimmed similar to TrueSeq:
   1. AdaptorRead1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
   2. AdaptorRead2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

These were the index primers, cDNA input, and PCR cycles used for each sample:

| Sample Type  | cDNA concentration | amount cDNA input (conc * 26 uL) | PCR cycles  for final step | Index_ID  | index  |
|--------------|--------------------|----------------------------------|------------|-----------|--------|
| #4 (Frag A)  | 0.0869             | 2.2594                           | 11         | E7500S-21 | GTTTCG |
| #5 (Frag A)  | 0.211              | 5.486                            | 9          | E7500S-14 | AGTTCC |
| #8 (Frag B)  | 0.0744             | 1.9344                           | 11         | E7500S-23 | GAGTGG |
| #9 (Frag B)  | 0.0232             | 0.6032                           | 11         | E7500S-18 | GTCCGC |
| #15 (Frag C) | 0.199              | 5.174                            | 9          | E7500S-25 | ACTGAT |
| #16 (Frag C) | 0.172              | 4.472                            | 9          | E7500S-20 | GTGGCC |
| #20 (Frag D) | 0.24               | 6.24                             | 9          | E7500S-13 | AGTCAA |
| #21 (Frag D) | 0.112              | 2.912                            | 9          | E7500S-15 | ATGTCA |
| #26 (Frag E) | 0.223              | 5.798                            | 9          | E7500S-22 | CGTACG |
| #27 (Frag E) | 0.188              | 4.888                            | 9          | E7500S-19 | GTGAAA |

Index primer full sequences from [oligo kit E7500S](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf?rev=2e735fd18b544d46b36ee0e88353ef5c&sc_lang=en-us&hash=CC77B45817715F3ED3A8F3B1953450EB):

| INDEX PRIMER | INDEX PRIMER SEQUENCE  | EXPECTED INDEX PRIMER SEQUENCE READ |
|---------------|------------------------|-----------------------------|
| NEBNext Index 13 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATTG_**TTGACT**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | AGTCAA | 
| NEBNext Index 14 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATAC_**GGAACT**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | AGTTCC | 
| NEBNext Index 15 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATTC_**TGACAT**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | ATGTCA | 
| NEBNext Index 16 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATGC_**GGACGG**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | CCGTCC | 
| NEBNext Index 18 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATGT_**GCGGAC**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | GTCCGC | 
| NEBNext Index 19 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATCG_**TTTCAC**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | GTGAAA | 
| NEBNext Index 20 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATAA_**GGCCAC**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | GTGGCC | 
| NEBNext Index 21 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATTC_**CGAAAC**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | GTTTCG | 
| NEBNext Index 22 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATTA_**CGTACG**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | CGTACG | 
| NEBNext Index 23 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATAT_**CCACTC**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | GAGTGG | 
| NEBNext Index 25 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATAT_**ATCAGT**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | ACTGAT | 
| NEBNext Index 27 Primer for Illumina (10 µM) | 5´-CAAGCAGAAGACGGCATACGAGATAA_**AGGAAT**_GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´ | ATTCCT | 