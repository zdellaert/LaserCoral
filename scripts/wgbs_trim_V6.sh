#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=200GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

# load modules needed

module load uri/main
module load fastp/0.23.2-GCC-11.2.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' LCM*R1*.fastq.gz))
R2_raw=($('ls' LCM*R2*.fastq.gz))

sample_name=$(basename $R1_raw | cut -d'_' -f1-2)

#make trimmed_V6_qc output folder
mkdir -p /home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V6_qc
mkdir -p /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V6/
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V6/"
qc_dir="/home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V6_qc"

for i in ${!R1_raw[@]}; do

sample_name=$(basename ${R1_raw[$i]} | cut -d'_' -f1-2)

  fastp \
  --in1 ${R1_raw[$i]} \
  --in2 ${R2_raw[$i]} \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --trim_poly_x \
  --thread 16 \
  --trim_front1 25 \
  --trim_front2 25 \
  --html ${qc_dir}/"$sample_name".fastp-trim.report.html \
  --json ${qc_dir}/"$sample_name".fastp-trim.report.json \
  --out1 ${output_dir}/"${R1_raw[$i]}".fastp-trim.fq.gz \
  --out2 ${output_dir}/"${R2_raw[$i]}".fastp-trim.fq.gz \
  2> ../scripts/outs_errs/"${stderr_PE_name}".fastp-trim.stderr

    echo "trimming of $sample_name: ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# go to directory with trimmed files

cd /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V6

# Create an array of fastq files to process
files=($('ls' *fastp-trim.fq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ${qc_dir}

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
