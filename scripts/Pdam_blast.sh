#!/bin/bash
#SBATCH --job-name="Pdam_blast"
#SBATCH -t 240:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH --mem=100GB
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --account=putnamlab
#SBATCH --nodes=2 --ntasks-per-node=24

module load BLAST+/2.15.0-gompi-2023a

cd ../references

makeblastdb -in Pdam_RefSeq_protein.faa -out blast_dbs/Pdam_prot -dbtype prot

blastp -query Pocillopora_acuta_HIv2.genes.pep.faa -db blast_dbs/Pdam_prot -out annotation/blastp_Pdam_out.tab -outfmt 6 -evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1

echo "Blast complete!" $(date)
