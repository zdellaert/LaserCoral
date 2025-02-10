#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=1 #split one task over multiple CPU
#SBATCH --mem=32GB
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load samtools/1.19.2

cd output_WGBS/dedup_V3

mkdir -p dedup_bam_subsets

for bamFile in *deduplicated.sorted.bam; do
        prefix=$(basename $bamFile .bam | cut -d'_' -f1-4)
        samtools index --threads 1 $bamFile
        samtools view -b $bamFile Pocillopora_acuta_HIv2___Sc0000000:1-1000000 > ${prefix}.deduplicated.sorted.sub.bam
        samtools index ${prefix}.deduplicated.sorted.sub.bam
done

mv *sub.bam* subsets/