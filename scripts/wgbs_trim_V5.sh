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
module load cutadapt/3.5-GCCcore-11.2.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' trimmed_V3_LCM*R1*.fastq.gz))
R2_raw=($('ls' trimmed_V3_LCM*R2*.fastq.gz))

#make trimmed_V5_qc output folder
mkdir -p /home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V5_qc
mkdir -p /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V5/
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V5/"
qc_dir="/home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V5_qc"

for i in ${!R1_raw[@]}; do
    cutadapt \
    -a "A{10}" \
    -A "A{10}" \
    -o  ${output_dir}trimmed_V5_${R1_raw[$i]} -p ${output_dir}trimmed_V5_${R2_raw[$i]} \
    ${R1_raw[$i]} ${R2_raw[$i]} \
    -q 20,20 --nextseq-trim=20 --minimum-length 40 --cores=20

    echo "trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# go to directory with trimmed files

cd /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V5

# Create an array of fastq files to process
files=($('ls' trimmed_V5*.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ${qc_dir}

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
