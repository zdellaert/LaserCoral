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

#load packages
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
