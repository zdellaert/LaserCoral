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
array_trim=($(ls trimmed_oligo*)) 

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
