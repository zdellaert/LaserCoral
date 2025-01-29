#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=200GB
#SBATCH -t 6:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

# load modules needed
module load uri/main
module load pigz/2.7-GCCcore-11.3.0
module load trimgalore/0.6.10

#make arrays of R1 and R2 reads trimmed above
to_trim=($('ls' trimmed_V3_LCM_3*.fastq.gz))

#make trimmed_V4_qc output folder
mkdir -p ../output_WGBS/trimmed_V4_qc/
mkdir -p /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V4

for i in ${!to_trim[@]}; do
    trim_galore ${to_trim[$i]} \
    --hardtrim5 100 \
    --cores 4 \
    --output_dir /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V4

    echo "trimming of ${to_trim[$i]} complete"
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# go to QC output directory
cd ../output_WGBS/trimmed_V4_qc/

# Create an array of fastq files to process
files=($('ls' /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V4/*100bp_5prime.fq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ./ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)

