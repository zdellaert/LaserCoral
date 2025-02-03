#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=128
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

# Make arrays of R1 and R2 reads
R1_raw=($(ls LCM_1_*R1*.fastq.gz))
R2_raw=($(ls LCM_1_*R2*.fastq.gz))

for i in ${!R1_raw[@]}; do

# base_name=${R1_raw[$i]%%_R1*.fastq.gz}

#     flexbar \
#     -r ${R1_raw[$i]} \
#     -p ${R2_raw[$i]} \
#     --adapter-preset TruSeq \
#     --adapter-min-overlap 3 \
#     --adapter-error-rate 0.1 \
#     --adapter-trim-end RIGHT \
#     --adapter-pair-overlap ON \
#     --qtrim TAIL --qtrim-threshold 20 \
#     --qtrim-format sanger \
#     --min-read-length 40 \
#     --htrim-right G \
#     --zip-output GZ \
#     --threads 120 \
#     -t ${base_name}_fbtrim

#     echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
# done

# for i in ${!R1_raw[@]}; do

base_name=${R1_raw[$i]%%_R1*.fastq.gz}

    flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    -a adapters.fasta \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end ANY \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \
    --zip-output GZ \
    --threads 120 \
    -t ${base_name}_fbtrim_test_ANY

    echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

for i in ${!R1_raw[@]}; do

base_name=${R1_raw[$i]%%_R1*.fastq.gz}

    flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    -a adapters.fasta \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end ANY \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \
    --htrim-left G \
    --zip-output GZ \
    --threads 120 \
    -t ${base_name}_fbtrim_test_RLG

    echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

# Create an array of fastq files to process
files=($('ls' *fbtrim_test*gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

multiqc *  #Compile MultiQC report from FastQC files 
