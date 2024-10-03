# LCM RNA-seq Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 10/2/2024

Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
RNA extractions: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Extractions/
RNA library prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Low-Input-RNA-Library-Prep/

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


