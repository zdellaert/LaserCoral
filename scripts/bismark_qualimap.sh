#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=1 #split one task over multiple CPU
#SBATCH --mem=32GB
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load qualimap/2.2.1

cd output_WGBS/dedup_V3

mkdir -p qualimap/bamqc

for bamFile in *deduplicated.sorted.bam; do
	prefix=$(basename $bamFile .bam)

	qualimap \
    --java-mem-size=29491M \
    bamqc \
     \
    -bam ${bamFile}  \
     \
    -p non-strand-specific \
    --collect-overlap-pairs \
    -outdir qualimap/bamqc/${prefix} \
    -nt 6
done
