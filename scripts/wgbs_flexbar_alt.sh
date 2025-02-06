#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=150GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

# Get the list of sample files and corresponding sample names
files=($(ls LCM*R1*.fastq.gz))
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.fastq.gz")
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/flexbar"
mkdir -p ${output_dir}

flexbar \
    -r ${sample_name}_R1_001.fastq.gz \
    -p ${sample_name}_R2_001.fastq.gz \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G --htrim-left G --htrim-min-length 5 \
    --max-uncalled 2 --htrim-error-rate 0.2 \
    --zip-output GZ \
    --threads 20 \
    -t "${output_dir}/${sample_name}_flexbar_alt"

# load modules needed
module load fastqc/0.12.1

#make trimmed_flexbar_qc output folder
mkdir ../output_WGBS/trimmed_flexbar_alt_qc/

# Create an array of fastq files to process
files=($('ls' ${output_dir}/${sample_name}_flexbar*gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/trimmed_flexbar_alt_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)
