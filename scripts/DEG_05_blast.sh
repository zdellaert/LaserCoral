#!/bin/bash
#SBATCH --job-name="DE_blast"
#SBATCH -t 500:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --mem=500GB
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --account=putnamlab
#SBATCH --nodes=2 --ntasks-per-node=24

module load BLAST+/2.15.0-gompi-2023a

cd ../output_RNA/differential_expression #set working directory
mkdir blast
cd blast

#nr database location andromeda: /data/shared/ncbi-db/.ncbirc 
# points to current location: cat /data/shared/ncbi-db/.ncbirc
# [BLAST]
# BLASTDB=/data/shared/ncbi-db/2024-11-10

blastp -query ../DEG_05_seqs.txt -db nr -out DEG_05_blast_results_restart.txt -outfmt 0 -evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 10

blastp -query ../DEG_05_seqs.txt -db nr -out DEG_05_blast_results_tab.txt -outfmt 6 -evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1

echo "Blast complete!" $(date)
