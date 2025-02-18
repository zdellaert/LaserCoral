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

multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

cat CpG.all.samps.5x_sorted.bed | awk '$4 ==10' > CpG.filt.all.samps.5x_sorted.bed
