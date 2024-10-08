#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=75GB
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
