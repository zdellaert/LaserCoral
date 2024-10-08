# LCM RNA-seq Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 10/5/2024

- Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
- RNA extractions: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Extractions/
- RNA library prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/

## Sequencing information

- Sequenced through Genohub Service Provider: Oklahoma Medical Research Foundation NGS Core
- Instrument: Illumina NovaSeq X Plus - 25B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Number of samples: 10
- Guaranteed number of pass filter PE reads/sample: 30M (15M in each direction)

## Make directory structure on Andromeda

```
mkdir /data/putnamlab/zdellaert/LaserCoral
cd /data/putnamlab/zdellaert/LaserCoral #Enter working directory

mkdir scripts #make folder for scripts
mkdir data_RNA  #make folder for raw data
mkdir output_RNA #make folder for output
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
array1=($('ls' *.fastq.gz)) 

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

### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/raw_qc), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/raw_qc_multiqc_report.html)

Okay! So there is definitely data! Yay!

And a lot of it!

And it isn't horrible quality! Yay!

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot1.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot2.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot3.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot4.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot5.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot6.png?raw=true">

But there is also a lot of duplication, overrepresentation of seqeunces, and adapter content. This isn't unexpected given the fact that this was low-input of degraded RNA extracted from fixed, LCM-d tissue. Let's see what we can do to interpret and filter from here.

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot_fastqc_adapter.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/raw_qc/multiqc_screenshots/screenshot_fastqc_overrep.png?raw=true">

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

## Trimming adapters and low-quality bases and short reads (< 20bp)

I am going to use [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and quality control

Example code with comments:

```
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \  # NEB AdaptorRead1 
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # NEB AdaptorRead2
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \ #output files
    input_R1.fastq.gz input_R2.fastq.gz \ #input files
    -q 20,20 \ #trims low-quality bases (score < 20) from the 3' end (first 20) and 5' (second 20) of the read
    --minimum-length 20 #after trimming, only keep a sequence if longer than 20 bp
```

```
nano scripts/cutadapt.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/data_RNA

# load modules needed
module load cutadapt/4.2-GCCcore-11.3.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' *R1*.fastq.gz))
R2_raw=($('ls' *R2*.fastq.gz))

for i in ${!R1_raw[@]}; do
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o trimmed_${R1_raw[$i]} -p trimmed_${R2_raw[$i]} \
    ${R1_raw[$i]} ${R2_raw[$i]} \
    -q 20,20 --minimum-length 20 --cores=20

    echo "trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# unload conflicting modules with modules needed below
module unload cutadapt/4.2-GCCcore-11.3.0
module unload GCCcore/11.3.0 Python/3.10.4-GCCcore-11.3.0 libffi/3.4.2-GCCcore-11.3.0

# load modules needed
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make trimmed_qc output folder
mkdir trimmed_qc/

# Make an array of fastq files to trim
array_trim=($('ls' trimmed*)) 

#run fastqc on trimmed data
for i in ${array_trim[@]}; do
    fastqc ${i} -o trimmed_qc/
    echo "fastqc ${i} done"
done

#Compile MultiQC report from FastQC files
multiqc trimmed_qc/  #Compile MultiQC report from FastQC files 

mv multiqc_report.html trimmed_qc/trimmed_qc_multiqc_report.html
mv multiqc_data trimmed_qc/trimmed_multiqc_data

echo "QC of trimmed data complete." $(date)
```
### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/trimmed_qc), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_qc/trimmed_qc_multiqc_report.html)

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_qc/multiqc_screenshots/screenshot1.png?raw=true">

Okay! So this trimming definitely removed the adapter sequences, and the QC data reflect that, but there are still oligo and primer sequences in the reads. 

## Trimming oligo and primer sequences from adapter-trimmed reads

1. Oligo Sequences
   1. NEBNext Template Switching Oligo: 5 ́-GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3 ́
   2. NEBNext Single Cell RT Primer Mix: 5 ́-AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTV-3 ́
   3. NEBNext Single Cell cDNA PCR Primer: 5 ́-AAGCAGTGGTATCAACGCAGAGT-3 ́

I am going to try to trim these using the following cutadapt parameters:

```
-a GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATGGG #trim Template Switching Oligo from R1
-a AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT #trim Single Cell RT Primer from R1
-a AAGCAGTGGTATCAACGCAGAGT #trim Single Cell cDNA PCR Primer from R1
-A GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATGGG #trim Template Switching Oligo from R2
-A AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT #trim  Single Cell RT Primer from R2
-A AAGCAGTGGTATCAACGCAGAGT #trim Single Cell cDNA PCR Primer from R2
``` 

```
nano scripts/cutadapt_oligo.sh #write script for further trimming pass and QC, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/data_RNA

# load modules needed
module load cutadapt/4.2-GCCcore-11.3.0

#make arrays of R1 and R2 reads
R1_trimmed=($('ls' trimmed_*R1*.fastq.gz))
R2_trimmed=($('ls' trimmed_*R2*.fastq.gz))

for i in ${!R1_trimmed[@]}; do
    cutadapt \
    -a GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATGGG \
    -a AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT \
    -a AAGCAGTGGTATCAACGCAGAGT \
    -A GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATGGG \
    -A AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTT \
    -A AAGCAGTGGTATCAACGCAGAGT \
    -o trimmed_oligo_${R1_trimmed[$i]} -p trimmed_oligo_${R2_trimmed[$i]} \
    ${R1_trimmed[$i]} ${R2_trimmed[$i]} \
    -q 20,20 --minimum-length 20 --cores=20

    echo "trimming of ${R1_trimmed[$i]} and ${R2_trimmed[$i]} complete"
done

# unload conflicting modules with modules needed below
module unload cutadapt/4.2-GCCcore-11.3.0
module unload GCCcore/11.3.0 Python/3.10.4-GCCcore-11.3.0 libffi/3.4.2-GCCcore-11.3.0

# load modules needed
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make trimmed_oligo_qc output folder
mkdir trimmed_oligo_qc/

# Make an array of fastq files to trim
array_trim=($('ls' trimmed_oligo*)) 

#run fastqc on trimmed_oligo data
for i in ${array_trim[@]}; do
    fastqc --threads 20 ${i} -o trimmed_oligo_qc/
    echo "fastqc ${i} done"
done

#Compile MultiQC report from FastQC files
multiqc trimmed_oligo_qc/  #Compile MultiQC report from FastQC files 

mv multiqc_report.html trimmed_oligo_qc/trimmed_oligo_qc_multiqc_report.html
mv multiqc_data trimmed_oligo_qc/trimmed_oligo_multiqc_data

echo "QC of trimmed_oligo data complete." $(date)
```

### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_RNA/trimmed_oligo_qc), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/trimmed_oligo_qc_multiqc_report.html)

Yay! I feel good about moving forward here to alignment. There is still duplication, but the non-biological overrepresented sequences appear to have been removed!

<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot1.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot2.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot3.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot4.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot5.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot6.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot7.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot8.png?raw=true">
<img width="800" alt="screenshot" src="https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/trimmed_oligo_qc/multiqc_screenshots/screenshot9.png?raw=true">

## Download Genome: [*Pocillopora acuta*](http://cyanophora.rutgers.edu/Pocillopora_acuta/)

Rutgers University Stephens et al. 2022 [Publication](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac098/6815755)

Obtain reference genome assembly and gff annotation file.

```
mkdir references/
cd references/ 

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.assembly.fasta.gz

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.genes.gff3.gz

gunzip Pocillopora_acuta_HIv2.assembly.fasta.gz #unzip genome file
gunzip Pocillopora_acuta_HIv2.genes.gff3.gz #unzip gff annotation file
```

### Convert GFF files to GTF format for Stringtie assembler

The gff3 file provided with the *Pocillopora acuta* genome is missing some features that are necessary for this pipeline (Stringtie specifically). I can correct the gff3 file and add those fields, but it is easier to convert the gff3 to a gtf file and automatically add those fields in the process. In order to do this, I am going to use the program [gffread](https://github.com/gpertea/gffread). Information and documentation about this package can be found on [the github examples page](https://github.com/gpertea/gffread/tree/master/examples).

```
nano scripts/gffread.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/references

# load modules needed
module load gffread/0.12.7-GCCcore-11.2.0

# "Clean" GFF file if necessary, then convert cleaned file into a GTF
# -E : remove any "non-transcript features and optional attributes"

gffread -E Pocillopora_acuta_HIv2.genes.gff3 -T -o Pocillopora_acuta_HIv2.gtf 

echo "Check how many fields are in each row of the file; currently there are rows with two different lenghts: 10 and 12"
awk '{print NF}' Pocillopora_acuta_HIv2.gtf | sort -u

# Use awk to add "gene_id = TRANSCRIPT_ID" to each of the rows that only have a transcript id listed (the non-transcript lines)
awk 'BEGIN {FS=OFS="\t"} {if ($9 ~ /transcript_id/ && $9 !~ /gene_id/) {match($9, /transcript_id "([^"]+)";/, a); $9 = $9 " gene_id \"" a[1] "\";"}; print}' Pocillopora_acuta_HIv2.gtf > Pocillopora_acuta_HIv2_modified.gtf

echo "Check how many fields are in each row of the file; Now all the rows should be the same length and only one value should be printed, 12"
awk '{print NF}' Pocillopora_acuta_HIv2_modified.gtf | sort -u

# remove the non-modified file
rm Pocillopora_acuta_HIv2.gtf

# rename the modified file
mv Pocillopora_acuta_HIv2_modified.gtf Pocillopora_acuta_HIv2.gtf
```

Script outputs:

- loaded 33730 genomic features from Pocillopora_acuta_HIv2.genes.gff3

- Check how many fields are in each row of the file; currently there are rows with two different lenghts: 10 and 12
    - 10
    - 12
- Check how many fields are in each row of the file; Now all the rows should be the same length and only one value should be printed, 12
    - 12

## HISAT2 Alignment

I will use [Hisat2](https://daehwankimlab.github.io/hisat2/manual/) to align the RNA-seq reads to the *P. acuta* genome

The libraries are paired and and UNstranded, since they were [prepared](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/) using a [template switching method](https://www.neb.com/en-us/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina) - not with the directional Ultra II NEB workflow.

See notes here: [strand-related settings for RNA-seq tools](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/)

```
hisat2 -p 16 \ #use 16 threads
    --time \ Print the wall-clock time required to load the index files and align the reads to stderr
    --dta \ #for input into Stringtie transcriptome assembly
    -q \ #fastq input files
    -x Pacuta_ref \ #index location 
    -1 ${read1} -2 ${read2} \ #input files, R1 and R2
    -S hisat2/${sample_name}.sam #output sam file
```

```
nano scripts/hisat2.sh #write script for alignment, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# load modules needed
module load HISAT2/2.2.1-gompi-2022a #Alignment to reference genome: HISAT2
module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

# index the reference genome, will write to a directory called Pacuta_ref
hisat2-build -f ../references/Pocillopora_acuta_HIv2.assembly.fasta ./Pacuta_ref

echo "Reference genome indexed. Starting alingment" $(date)

# make the output directory if it does not exist (-p checks for this)
mkdir -p hisat2

# call the oligo-trimmed sequences into an array
array=(../data_RNA/trimmed_oligo*_R1_001.fastq.gz)

# align the files to the indexed genome using hisat2
for read1 in ${array[@]}; do
    
    # extract the sample name of R1 files as LCM_##
    sample_name=$(basename $read1 | sed 's/.*trimmed_\([A-Za-z0-9]*_[0-9]*\).*/\1/')
    
    # Define the corresponding reverse read file (R2)
    read2=${read1/_R1_/_R2_}
    
    # perform alignment
    hisat2 -p 16 --time --dta -q -x Pacuta_ref -1 ${read1} -2 ${read2} -S hisat2/${sample_name}.sam
    echo "${sample_name} aligned!"

    # sort the sam file into a bam file
    samtools sort -@ 8 -o hisat2/${sample_name}.bam hisat2/${sample_name}.sam
    echo "${sample_name} bam-ified!"
    
    # index bam file , creating a .bai file which is nice for viewing in IGB
    samtools index hisat2/${sample_name}.bam hisat2/${sample_name}.bai
    
    # remove sam file to save disk space
    rm hisat2/${sample_name}.sam
done

# move the reference index files into the hisat2 directory
mv Pacuta_ref.* hisat2/

#  Calculate mapping percentages
for i in hisat2/*.bam; do
    echo "${i}" >> mapped_reads_counts_Pacuta
    samtools flagstat ${i} | grep "mapped (" >> hisat2/mapped_reads_counts_Pacuta
done
```

Average mapping rate: **81.71%**
- min 79.72%
- max 83.64%
- average 81.71%

Average primary mapping rate: **70.55%** (The percentage of reads mapped in their primary alignment (without secondary alignments))
- min 65.70%
- max 73.58%
- average 70.55%

### Stranded HISAT2 Alignment

I will use [Hisat2](https://daehwankimlab.github.io/hisat2/manual/) to align the RNA-seq reads to the *P. acuta* genome

The libraries are paired and could porentially be stranded, since they were [prepared](https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/) using a [template switching method](https://www.neb.com/en-us/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina).

See notes here: [strand-related settings for RNA-seq tools](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/)

```
hisat2 -p 16 \ #use 16 threads
    --time \ Print the wall-clock time required to load the index files and align the reads to stderr
    --dta \ #for input into Stringtie transcriptome assembly
    -q \ #fastq input files
    **--rna-strandness RF** \
    -x Pacuta_ref \ #index location 
    -1 ${read1} -2 ${read2} \ #input files, R1 and R2
    -S hisat2/${sample_name}.sam #output sam file
```

```
nano scripts/stranded_hisat2.sh #write script for alignment, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# load modules needed
module load HISAT2/2.2.1-gompi-2022a #Alignment to reference genome: HISAT2
module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

# index the reference genome, will write to a directory called Pacuta_ref
hisat2-build -f ../references/Pocillopora_acuta_HIv2.assembly.fasta ./Pacuta_ref

echo "Reference genome indexed. Starting alingment" $(date)

# make the output directory if it does not exist (-p checks for this)
mkdir -p stranded_hisat2

# call the oligo-trimmed sequences into an array
array=(../data_RNA/trimmed_oligo*_R1_001.fastq.gz)

# align the files to the indexed genome using hisat2
for read1 in ${array[@]}; do
    
    # extract the sample name of R1 files as LCM_##
    sample_name=$(basename $read1 | sed 's/.*trimmed_\([A-Za-z0-9]*_[0-9]*\).*/\1/')
    
    # Define the corresponding reverse read file (R2)
    read2=${read1/_R1_/_R2_}
    
    # perform alignment
    hisat2 -p 16 --time --dta -q --rna-strandness RF -x Pacuta_ref -1 ${read1} -2 ${read2} -S stranded_hisat2/${sample_name}.sam
    echo "${sample_name} aligned!"

    # sort the sam file into a bam file
    samtools sort -@ 8 -o stranded_hisat2/${sample_name}.bam stranded_hisat2/${sample_name}.sam
    echo "${sample_name} bam-ified!"
    
    # index bam file , creating a .bai file which is nice for viewing in IGB
    samtools index stranded_hisat2/${sample_name}.bam stranded_hisat2/${sample_name}.bai
    
    # remove sam file to save disk space
    rm stranded_hisat2/${sample_name}.sam
done

# move the reference index files into the stranded_hisat2 directory
mv Pacuta_ref.* stranded_hisat2/

#  Calculate mapping percentages
for i in stranded_hisat2/*.bam; do
    echo "${i}" >> stranded_hisat2/mapped_reads_counts_Pacuta
    samtools flagstat ${i} | grep "mapped (" >> stranded_hisat2/mapped_reads_counts_Pacuta
done
```

### Interpretation of unstranded vs. stranded hisat2 run

- The output and error files are the exact same , so the alignments did not seem to be affected by the stradedness.
- Upon inspecting the bam files to see if the information about alignments was affected, there do not appear to be any differences between the stranded vs. unstranded runs.

#### RSeQC: infer_experiment.py - infer strandedness

I found [this tool](https://rseqc.sourceforge.net/#infer-experiment-py) to analyse aligned bam files to examine for stradedness information. Running it on both the stranded and unstrandedd hisat2 alignments:

```
cd strandedness
pip install RSeQC
cp stranded_test_trimmed_oligo_trimmed_LCM_15_S40_R1_001/Pocillopora_acuta_HIv2.bed  . #copy bed file made using other tool

infer_experiment.py -r Pocillopora_acuta_HIv2.bed -i ../output_RNA/hisat2/LCM_15.bam
```

```
Reading reference gene model Pocillopora_acuta_HIv2.bed ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.0002
Fraction of reads explained by "1++,1--,2+-,2-+": 0.5102
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4897
```
```
infer_experiment.py -r Pocillopora_acuta_HIv2.bed -i ../output_RNA/stranded_hisat2/LCM_15.bam
```

```
Reading reference gene model Pocillopora_acuta_HIv2.bed ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.0002
Fraction of reads explained by "1++,1--,2+-,2-+": 0.5102
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4897
```

For both the bam files, the results were the same. The reads are distributed almost equally across both orientations. This indicates that the library prep method I used is **unstranded**, even though it was using a directional selection (3' poly-A primer) method. This is what I thought, but wanted to double check.

## Assembly with Stringtie

I will use [Stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual) to perform reference-guided assembly of the RNA-seq data. I am using the simplified stringtie protocol with the "stringtie -eB" option:

<img width="800" alt="stringtie_workflow" src="https://github.com/zdellaert/LaserCoral/blob/main/code/images/stringtie_workflow.png?raw=true">

Note: Use of -e (Expression estimation mode) means stringtie will only estimate the abundance of given reference transcripts (only genes from the genome, no novel genes). If you are interested in novel genes or splice variants, see Stringtie documentation:

> StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file. With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored.

```
nano scripts/stringtie.sh #make script for assembly, enter text in next code chunk
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# move into stringtie directory
module load StringTie/2.1.4-GCC-9.3.0

# make the output directory if it does not exist (-p checks for this)
mkdir -p stringtie

# call the hisat2 bam files into an array
array=(hisat2/*.bam)

for i in "${array[@]}"; do 
    sample_name=$(echo "$i" | awk -F'[/.]' '{print $2}')

    # -p 16 : use 16 cores
    # -e : exclude novel genes
    # -B : create Ballgown input files for downstream analysis
    # -v : enable verbose mode
    # -G : gtf annotation file
    # -A : output name for gene abundance estimate files
    # -o : output name for gtf file

    stringtie -p 16 -e -B -v \
        -G ../references/Pocillopora_acuta_HIv2.gtf \
        -A stringtie/"${sample_name}".gene_abund.tab \
        -o stringtie/"${sample_name}".gtf \
        "$i" #input bam file

    echo "StringTie assembly for seq file ${i}" $(date)
done

echo "StringTie assembly COMPLETE" $(date)
```

## Generate gene count matrix using [prepDE.py script from Stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

Download script from [stringtie website](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3) or [github repository](https://github.com/gpertea/stringtie/blob/master/prepDE.py3). I am using the python3 version, but this and the original version (prepDE.py) are very similar and should give the exact same result. I am using [this input file format](https://ccb.jhu.edu/software/stringtie/dl/sample_lst.txt).

```
cd scripts
wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
nano prepDE.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# load modules needed
module load StringTie/2.1.4-GCC-9.3.0 #stringtie module, includes Python 3.8.5

# move into stringtie directory
cd stringtie

# make input file
for filename in *.gtf; do
    sample_name=$(echo "$filename" | awk -F'[.]' '{print $1}')

    echo $sample_name $PWD/$filename
done > listGTF.txt

#Compile the gene count matrix
python ../../scripts/prepDE.py3 -g LCM_RNA_gene_count_matrix.csv -i listGTF.txt

echo "Gene count matrix compiled." $(date)
```

Woohoo! [Gene count matrix complete.](https://github.com/zdellaert/LaserCoral/blob/main/output_RNA/stringtie/LCM_RNA_gene_count_matrix.csv)

## GeneExt to Extend Genes to include putative 3' UTRs

[GeneExt](https://github.com/sebepedroslab/GeneExt/tree/main) is a program from the Sebé-Pedrós Lab which has a great manual that explains this probem with 3' biased sequencing approaches in non-model organism in great detail and with helpful nuance. The program allows the user to create a modified GTF/GFF file based on a BAM file of mapping data. It appends 3' UTRs to the genes in a way that is informed by the mapping of reads, so the length of the 3' UTRs added can change dynamically for each gene based on where reads are acutally mapping.

See manual [here](https://github.com/sebepedroslab/GeneExt/blob/main/Manual.md)

### Installation

I have done this on andromeda, in interactive mode. 

```
interactive -c 24

cd /data/putnamlab/zdellaert/
module load Miniconda3/4.9.2
cd GeneExt

# create environment
conda env create -n geneext -f environment.yaml

# update environment
conda env update -n geneext -f environment.yaml
```

### Combine all bam files

Combine all 10 bam files together

```
nano scripts/bammerge.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

#load packages
module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

#use samtools merge to merge all the files
samtools merge hisat2/merge.bam hisat2/*.bam
```

### Running GeneExt

```
nano scripts/GeneExt.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/references #set working directory

# Use the gtf file run in Stringtie above:

# the gtf has to have the gene and transcript IDs be different, so append -T to each transcript_id and save to a new file
sed 's/transcript_id "\([^"]*\)"/transcript_id "\1-T"/g' Pocillopora_acuta_HIv2.gtf > Pocillopora_acuta_HIv2_modified.gtf

# activate environment
source ~/.bashrc
module load Miniconda3/4.9.2
conda activate geneext

mkdir -p geneext
cd geneext

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand

python /data/putnamlab/zdellaert/GeneExt/geneext.py --verbose=3 \
    -g ../Pocillopora_acuta_HIv2_modified.gtf \
    -b ../../output_RNA/hisat2/merge.bam \
    -o Pocillopora_acuta_GeneExt.gtf \
    -j 18 --clip_strand both
```


## Downstream analysis thoughts

1. What level of filtering do we want pre-DEseq?
2. Any reason to do WGCNA or glmmseq instead of DEseq?
   1. glmmseq would allow us to incorporate possible batch effects of LCM date
3. Want to look into
   1. Differential expression of biomineralization toolkit
   2. Differential expression of transcription factors
   3. Differential expression of single cell marker genes
