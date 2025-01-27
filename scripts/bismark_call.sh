#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=50 #split one task over multiple CPU
#SBATCH --mem=500GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load parallel/20240822
module load uri/main
module load Bismark/0.23.1-foss-2021b

cd output_WGBS/dedup_V3

find *deduplicated.bismark.cov.gz | \
xargs -n 1 basename -s _pe.deduplicated.bismark.cov.gz | \
parallel -j 48 coverage2cytosine \
--genome_folder ../../references/ \
-o {} \
--merge_CpG \
--zero_based \
{}_pe.deduplicated.bismark.cov.gz
