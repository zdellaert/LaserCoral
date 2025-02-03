#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --mem=50GB
#SBATCH -t 2:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%A_%a"
#SBATCH --output=../scripts/outs_errs/"%x_output.%A_%a"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
#SBATCH --array=1-10  # Array job for 10 samples

# Define input and output directories
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V7/"
qc_dir="/home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V7_qc"

# Create necessary directories
mkdir -p ${output_dir}
mkdir -p ${qc_dir}

# Get list of input files
R1_raw=($('ls' LCM*R1*.fastq.gz))
R2_raw=($('ls' LCM*R2*.fastq.gz))

# Get the current sample index based on SLURM_ARRAY_TASK_ID (1-based index)
i=$((SLURM_ARRAY_TASK_ID - 1))

# Extract sample name
sample_name=$(basename ${R1_raw[$i]} | cut -d'_' -f1-2)

# Run flexbar for trimming
flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \
    --zip-output GZ \
    --threads 40 \
    -t ${output_dir}/"${sample_name}"_fbtrim

echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"

# Load modules
module load fastqc/0.12.1

# Run FastQC on trimmed output
fastqc ${output_dir}/${sample_name}_fbtrim_1.fastq.gz -o ${qc_dir}
fastqc ${output_dir}/${sample_name}_fbtrim_2.fastq.gz -o ${qc_dir}

echo "FastQC of ${sample_name} complete"

# MultiQC only runs once for the whole dataset, so it should not be part of the array job
