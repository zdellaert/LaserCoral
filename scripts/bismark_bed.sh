#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

cd output_WGBS/dedup_V3

# Bedgraphs for 5X coverage 

for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${sample}"_5x_sorted.bedgraph
done

# Bedgraphs for 10X coverage 

for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${sample}"_10x_sorted.bedgraph
done
