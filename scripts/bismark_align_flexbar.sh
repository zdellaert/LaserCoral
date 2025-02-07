#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=24 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

# Set directories and files
reads_dir="data_WGBS/" #directory containing trimmed data to align
genome_folder="references/" #directory containing original unmodified genome fasta and bismark Bisulfite_Genome directory

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_flexbar"
mkdir -p $output_dir

# Get the list of sample files and corresponding sample names
files=(${reads_dir}LCM_*_flexbar_1.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_flexbar_1.fastq.gz")

# Define log files for stdout and stderr
stdout_log="${output_dir}${sample_name}_stdout.log"
stderr_log="${output_dir}${sample_name}_stderr.log"

# Run Bismark align
bismark \
    -genome ${genome_folder} \
    -p 4 \
    -score_min L,0,-1.0 \
    -1 ${reads_dir}${sample_name}_flexbar_1.fastq.gz \
    -2 ${reads_dir}${sample_name}_flexbar_2.fastq.gz \
    -o ${output_dir} \
    --temp_dir ${output_dir} \
    --basename ${sample_name} \
    2> "${output_dir}/${sample_name}-bismark_summary.txt"

# Define directories
summary_file="${output_dir}/parameter_comparison_summary.csv"

# Initialize summary file
echo "Sample,Score_Min,Alignment_Rate" > ${summary_file}

for file in ${output_dir}/*_report.txt; do
    # Extract sample name and from directory name
    sample_name=$(basename "$file" | cut -d'_' -f1-3)
    score_min="L0-1.0"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
done
