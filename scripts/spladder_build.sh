#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/spladder

spladder_dir="output_RNA/splicing"
out_dir="/scratch3/workspace/zdellaert_uri_edu-shared/spladder_out"

mkdir -p ${out_dir}

spladder build --bams ${spladder_dir}/alignments.txt \
               --annotation ${spladder_dir}/Pocillopora_acuta_HIv2.gtf \
               --outdir ${out_dir}
