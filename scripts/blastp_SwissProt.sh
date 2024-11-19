#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral #set working directory

module load blast-plus/2.14.1

cd references/
mkdir annotation

fasta="Pocillopora_acuta_HIv2.genes.pep.faa"

blastp \
-query $fasta \
-db blast_dbs/uniprot_sprot_r2024_10_02 \
-out annotation/blastp_SwissProt_out.tab \
-evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6

echo "Blast complete!" $(date)
