## Hisat2-Bismark Alignment of V3 Trimmed Reads

**I had to delete the bowtie Bisulfite_Genome folder because the command wouldn't make a hisat version, it just kept copying the original one**

### Prepare Genome

```
nano scripts/bismark_genome_hisat.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 2:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load HISAT2/2.2.1-gompi-2022a

### Prepare Genome

mkdir -p references/Bisulfite_Genome_hisat2
cd references/Bisulfite_Genome_hisat2
cp ../Pocillopora_acuta_HIv2.assembly.fasta .

bismark_genome_preparation --hisat2 --verbose --parallel 10 ./
```

### align

```
nano scripts/bismark_align_hisat_V3.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=150GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load HISAT2/2.2.1-gompi-2021b

# Set directories and files
reads_dir="data_WGBS/" #directory containing trimmed data to align
genome_folder="references/Bisulfite_Genome_hisat2/" #directory containing original unmodified genome fasta and bismark Bisulfite_Genome directory

mkdir -p /scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_hisat_V3 #make output directory if it does not exist

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_hisat_V3"
checkpoint_file="/scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_hisat_V3/completed_samples.log"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}trimmed_V3_LCM_*R1_001.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.fastq.gz")

# Check if the sample has already been processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Define log files for stdout and stderr
stdout_log="${output_dir}${sample_name}_stdout.log"
stderr_log="${output_dir}${sample_name}_stderr.log"

# Run Bismark align
bismark \
    --hisat2 \
    -genome ${genome_folder} \
    -p 4 \
    -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
    -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
    -o ${output_dir} \
    --temp_dir ${output_dir} \
    --basename ${sample_name} \
    2> "${output_dir}/${sample_name}-bismark_summary.txt"

# Check if the command was successful
if [ $? -eq 0 ]; then
    # Append the sample name to the checkpoint file
    echo ${sample_name} >> ${checkpoint_file}
    echo "Sample ${sample_name} processed successfully."
else
    echo "Sample ${sample_name} failed. Check ${stderr_log} for details."
fi

# Define directories
summary_file="${output_dir}/parameter_comparison_summary.csv"

# Initialize summary file
echo "Sample,Score_Min,Alignment_Rate" > ${summary_file}

for file in ${output_dir}/*_report.txt; do
    # Extract sample name and from directory name
    sample_name=$(basename "$file" | cut -d'_' -f1-5)
    short_name=$(basename "$file" | cut -d'_' -f3-4)
    score_min="hisat"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${short_name},${score_min},${mapping}" >> ${summary_file}
done
```

### hisat interpretation

Hisat with default parameters had decreased alignment rates. Could try making the hisat parameters less stringent but there is less documentation on this. 

<img src="../08-Bismark-Alignment-Assesment-images/alignment_bismark_bowtie_vs_hisat.png?raw=true" height="400">


## Pdam genome test

### Download pdam genome

```
mkdir references_Pdam/ 
cd references_Pdam/ 

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.fna.gz #download reference genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/704/095/GCF_003704095.1_ASM370409v1/GCF_003704095.1_ASM370409v1_genomic.gff.gz #download genome annotation file

gunzip GCF_003704095.1_ASM370409v1_genomic.fna.gz #unzip genome file
gunzip GCF_003704095.1_ASM370409v1_genomic.gff.gz #unzip gff annotation file

mv GCF_003704095.1_ASM370409v1_genomic.fna GCF_003704095.1_ASM370409v1_genomic.fasta #rename to .fasta 
```

### Prepare Genome

See: https://marineomics.github.io/FUN_02_DNA_methylation.html and https://felixkrueger.github.io/Bismark/bismark/genome_preparation/

See RNA seq bioniformatics markdown for genome download info.

```
nano scripts/bismark_genome_Pdam.sh 
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

cd references_Pdam/

bismark_genome_preparation --verbose --parallel 10 ./
```

### Bismark Alignment of Flexbar Trimmed Reads to Pdam Genome

```
nano scripts/bismark_align_flexbar_Pdam.sh 
```

```
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
genome_folder="references_Pdam/" #directory containing original unmodified genome fasta and bismark Bisulfite_Genome directory

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_flexbar_Pdam"
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
```

### interpretation:

Worse than the alignment with same parameters to P. acuta, so sticking with P. acuta. But a good test!

## Bismark Alignment of Flexbar Trimmed Reads

```
nano scripts/bismark_align_flexbar.sh 
```

```
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
```

## Deduplication and methylation extraction

```
nano scripts/bismark_flexbar_dedup_call.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load parallel/20240822
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/output_WGBS/align_flexbar"
ref_dir="/project/pi_hputnam_uri_edu/zdellaert/LaserCoral/references"

find ${output_dir}/*.bam | \
xargs -n 1 basename -s .bam | \
parallel -j 10 deduplicate_bismark \
--bam \
--paired \
--output_dir ${output_dir} \
${output_dir}/{}.bam

cd ${output_dir}

bismark_methylation_extractor \
--bedGraph \
--counts \
--comprehensive \
--merge_non_CpG \
--multicore 28 \
--buffer_size 75% \
*deduplicated.bam

module load all/MultiQC/1.12-foss-2021b

bismark2report
bismark2summary *pe.bam

bam2nuc --genome_folder ${ref_dir} *_pe.deduplicated.bam

multiqc .
```

There was very little difference between the final results of cutadapt vs. flexbar-trimmed reads using bismark.

## BWA alignment of flexbar-trimmed reads

Ran methylseq-bwameth on the flexbar trimmed reads, and determined that the alignment was less successful than with the cutadapt-trimmed reads. I will continue with the cutadapt trimming as my main trimming.

```
nano data_WGBS/LCM_methylseq_input_fb.csv
```

```
sample,fastq_1,fastq_2,genome
LCM_1,data_WGBS/LCM_1_S1_flexbar_1.fastq.gz,data_WGBS/LCM_1_S1_flexbar_2.fastq.gz,
LCM_3,data_WGBS/LCM_3_S2_flexbar_1.fastq.gz,data_WGBS/LCM_3_S2_flexbar_2.fastq.gz,
LCM_11,data_WGBS/LCM_11_S3_flexbar_1.fastq.gz,data_WGBS/LCM_11_S3_flexbar_2.fastq.gz,
LCM_12,data_WGBS/LCM_12_S4_flexbar_1.fastq.gz,data_WGBS/LCM_12_S4_flexbar_2.fastq.gz,
LCM_17,data_WGBS/LCM_17_S7_flexbar_1.fastq.gz,data_WGBS/LCM_17_S7_flexbar_2.fastq.gz,
LCM_18,data_WGBS/LCM_18_S8_flexbar_1.fastq.gz,data_WGBS/LCM_18_S8_flexbar_2.fastq.gz,
LCM_24,data_WGBS/LCM_24_S5_flexbar_1.fastq.gz,data_WGBS/LCM_24_S5_flexbar_2.fastq.gz,
LCM_25,data_WGBS/LCM_25_S6_flexbar_1.fastq.gz,data_WGBS/LCM_25_S6_flexbar_2.fastq.gz,
LCM_32,data_WGBS/LCM_32_S9_flexbar_1.fastq.gz,data_WGBS/LCM_32_S9_flexbar_2.fastq.gz,
LCM_33,data_WGBS/LCM_33_S10_flexbar_1.fastq.gz,data_WGBS/LCM_33_S10_flexbar_2.fastq.gz,
```

```
nano scripts/methylseq_fb.sh 
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=12
#SBATCH --mem=200GB
#SBATCH -t 4:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

## Load Nextflow and Apptainer environment modules
module purge
module load nextflow/24.10.3
module load apptainer/latest

## Set Nextflow directories to use scratch
export NXF_WORK=/scratch3/workspace/zdellaert_uri_edu-shared/nextflow_work
export NXF_TEMP=/scratch3/workspace/zdellaert_uri_edu-shared/nextflow_temp
export NXF_LAUNCHER=/scratch3/workspace/zdellaert_uri_edu-shared/nextflow_launcher

 nextflow run nf-core/methylseq \
   --input ./data_WGBS/LCM_methylseq_input_fb.csv \
   --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_bwa_test \
   --aligner bwameth \
   --methyl_kit \
   --igenomes_ignore \
   --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
   --relax_mismatches \
   --save_align_intermeds \
   --skip_trimming \
   --email zdellaert@uri.edu \
   -profile unity \
   -name methylseq_bwa_test
```