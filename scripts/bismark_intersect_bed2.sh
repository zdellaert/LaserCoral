#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

cd output_WGBS/dedup_V3

for i in *5x_sorted.bedgraph
do
   sample=$(basename ${i} _5x_sorted.bedgraph)
   bedtools map -a ../../references/Pocillopora_acuta_HIv2.gtf.cleaned.bed -b  ${i} -c 4 -o mean >  ${sample}_gene_body_methylation.txt

done
