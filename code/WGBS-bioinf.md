# LCM RNA-seq Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 1/16/2025

- Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
-  extractions: 
-  library prep: 

## Sequencing information

- Sequenced through Genohub Service Provider: Oklahoma Medical Research Foundation NGS Core

## Make directory structure on Unity within LaserCoral repo folder

```
cd /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/
mkdir data_WGBS
cd data_WGBS #Enter working directory
```

## Check integrity of data transfer

raw data is located in `/work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/`

```
md5sum /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*gz > URI_WGBS.md5
cat /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*md5 > Genohub_WGBS.md5

#use diff command to see if there is a differnece between the checksums
# using -w to ignore spaces (the URI file is formatted slightly differently) 

diff -w Genohub_WGBS.md5 URI_WGBS.md5

# no output, so no difference between the md5s. verified by inspecting manually.
```

Data appears to have been transferred successfully from genohub.

## Symlink raw data files into data_RNA

```
ln -s /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*gz /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
```

## QC raw files

```
nano ../scripts/wgbs_raw_qc.sh #write script for QC, enter text in next code chunk
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS


# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

#make raw_qc output folder
mkdir ../output_WGBS/raw_qc_WGBS/


# Create an array of fastq files to process
files=($('ls' *.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/raw_qc_WGBS/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

#Compile MultiQC report from FastQC files
multiqc ../output_WGBS/raw_qc_WGBS/  

mv multiqc_report.html ../output_WGBS/raw_qc_WGBS/raw_qc_multiqc_report.html
mv multiqc_data ../output_WGBS/raw_qc_WGBS/raw_multiqc_data

echo "Initial QC of Seq data complete." $(date)
```

### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_WGBS/raw_qc_WGBS), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/main/output_WGBS/raw_qc_WGBS/raw_qc_multiqc_report.html)

Okay! So there is definitely data! Yay!

And a lot of it!

And it isn't horrible quality! Yay!

But there is a lot of adapter content, duplication, and the GC content is wacky.

## Trimming adapters and low-quality bases and short reads (< 20bp)

From the Zymo trio protocol:

"Libraries should be trimmed to remove any adapter sequence. No other trimming is required. Use the following sequences to trim the adapters:"

Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

I am going to use [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and quality control

Example code with comments:

```
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # Zymo Adapter Read1 
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # Zymo Adapter Read2
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