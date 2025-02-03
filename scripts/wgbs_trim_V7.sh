#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS


# Define input and output directories
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V7/"
qc_dir="/home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V7_qc"

# Create necessary directories
mkdir -p ${output_dir}
mkdir -p ${qc_dir}

#make arrays of R1 and R2 reads
R1_raw=($('ls' LCM*R1*.fastq.gz))
R2_raw=($('ls' LCM*R2*.fastq.gz))


for i in ${!R1_raw[@]}; do

sample_name=$(basename ${R1_raw[$i]} | cut -d'_' -f1-2)

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
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# go to directory with trimmed files

cd ${output_dir}

# Create an array of fastq files to process
files=($('ls' *.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ${qc_dir}

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
