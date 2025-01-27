#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b

cd output_WGBS/dedup_V3

bismark_methylation_extractor \
--bedGraph \
--counts \
--comprehensive \
--merge_non_CpG \
--multicore 28 \
--buffer_size 75% \
*deduplicated.bam
