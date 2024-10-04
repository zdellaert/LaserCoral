#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error=../"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
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
