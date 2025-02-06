#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --mem=200GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

module load bedtools2/2.31.1

mkdir -p output_WGBS/dedup_V3/bam2fq
cd output_WGBS/dedup_V3/bam2fq

bedtools bamtofastq -i ../trimmed_V3_LCM_24_S5_pe.deduplicated.bam -fq dedup_V3_LCM_24_R1_bam2fq.fq -fq2 dedup_V3_LCM_24_R2_bam2fq.fq
bedtools bamtofastq -i ../trimmed_V3_LCM_33_S10_pe.deduplicated.bam -fq dedup_V3_LCM_33_R1_bam2fq.fq -fq2 dedup_V3_LCM_33_R2_bam2fq.fq

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1

# go to directory with trimmed files

# Create an array of fastq files to process
files=($('ls' *bam2fq*)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o . && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)
