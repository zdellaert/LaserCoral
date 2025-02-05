# LCM RNA-seq Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 1/16/2025

- Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
-  extractions: 
-  library prep: 

## Sequencing information

- Sequenced through Genohub Service Provider: Oklahoma Medical Research Foundation NGS Core

## Make directory structure on Unity within LaserCoral repo folder

```
cd /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/
mkdir data_WGBS
cd data_WGBS #Enter working directory
```

## Check integrity of data transfer

raw data is located in `/work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/`

```
md5sum /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*gz > URI_WGBS.md5
cat /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*md5 > Genohub_WGBS.md5

#use diff command to see if there is a differnece between the checksums
# using -w to ignore spaces (the URI file is formatted slightly differently) 

diff -w Genohub_WGBS.md5 URI_WGBS.md5

# no output, so no difference between the md5s. verified by inspecting manually.
```

Data appears to have been transferred successfully from genohub.

## Symlink raw data files into data_RNA

```
ln -s /work/pi_hputnam_uri_edu/20250113_LaserCoraL_WGBS/*gz /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
```

## QC raw files

```
nano ../scripts/wgbs_raw_qc.sh #write script for QC, enter text in next code chunk
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS


# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

#make raw_qc output folder
mkdir ../output_WGBS/raw_qc_WGBS/


# Create an array of fastq files to process
files=($('ls' *.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/raw_qc_WGBS/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

#Compile MultiQC report from FastQC files
multiqc ../output_WGBS/raw_qc_WGBS/  

mv multiqc_report.html ../output_WGBS/raw_qc_WGBS/raw_qc_multiqc_report.html
mv multiqc_data ../output_WGBS/raw_qc_WGBS/raw_multiqc_data

echo "Initial QC of Seq data complete." $(date)
```

### Interpretation of QC data

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_WGBS/raw_qc_WGBS), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/main/output_WGBS/raw_qc_WGBS/raw_qc_multiqc_report.html)

Okay! So there is definitely data! Yay!

And a lot of it!

And it isn't horrible quality! Yay!

But there is a lot of adapter content, duplication, and the GC content is wacky.

## Trimming adapters and low-quality bases and short reads (< 20bp)

From the Zymo trio protocol:

"Libraries should be trimmed to remove any adapter sequence. No other trimming is required. Use the following sequences to trim the adapters:"

Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

I am going to use [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and quality control

Example code with comments:

```
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # Zymo Adapter Read1 
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # Zymo Adapter Read2
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \ #output files
    input_R1.fastq.gz input_R2.fastq.gz \ #input files
    -q 20,20 \ #trims low-quality bases (score < 20) from the 3' end (first 20) and 5' (second 20) of the read
    --minimum-length 20 #after trimming, only keep a sequence if longer than 20 bp
```

```
nano scripts/wgbs_cutadapt.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS


# load modules needed

module load uri/main
module load cutadapt/3.5-GCCcore-11.2.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' *R1*.fastq.gz))
R2_raw=($('ls' *R2*.fastq.gz))

for i in ${!R1_raw[@]}; do
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o trimmed_${R1_raw[$i]} -p trimmed_${R2_raw[$i]} \
    ${R1_raw[$i]} ${R2_raw[$i]} \
    -q 20,20 --minimum-length 20 --cores=20

    echo "trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# unload conflicting modules with modules needed below
module unload cutadapt/3.5-GCCcore-11.2.0

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

#make trimmed_qc output folder
mkdir ../output_WGBS/trimmed_qc/

# Create an array of fastq files to process
files=($('ls' trimmed*.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/trimmed_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ../output_WGBS/trimmed_qc/

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
```

### Trimming done, QC looks good!

Will add screenshots soon.

## Bismark Alignment to P. acuta genome

https://felixkrueger.github.io/Bismark/

So, the version of nextcore that methylseq wants is >24.10, and unity only has 24.04. Going to run Bismark without the methylseq pipeline.

### Prepare Genome

See: https://marineomics.github.io/FUN_02_DNA_methylation.html and https://felixkrueger.github.io/Bismark/bismark/genome_preparation/

See RNA seq bioniformatics markdown for genome download info.

```
nano scripts/bismark_genome.sh 
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

cd references

bismark_genome_preparation --verbose --parallel 10 ./
```

### Aligning Reads

Trying an array job based on https://sr320.github.io/tumbling-oysters/posts/33-bismark-array/

I did use different alignment parameters from the above link but used the scripting format. For details about the bismark options see here: https://github.com/FelixKrueger/Bismark/blob/f88314914725242c59289a534271b66cf9d461e3/docs/options/alignment.md

Importantly: I have directional libraries, so I used the default directional setting

#### Testing score min parameters

based on https://sr320.github.io/tumbling-oysters/posts/33-bismark-array/

```
nano scripts/bismark_align_paramtest.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=100GB
#SBATCH -t 00:30:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

# Set directories and files
reads_dir="data_WGBS/"
genome_folder="references/"

mkdir -p output_WGBS/align_paramtest #make output directory if it does not exist

output_dir="output_WGBS/align_paramtest"
checkpoint_file="output_WGBS/align_paramtest/completed_samples.log"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}trimmed_LCM_*R1_001.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.fastq.gz")

# Check if the sample has already been processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Define log files for stdout and stderr
stdout_log="${output_dir}/${sample_name}_stdout.log"
stderr_log="${output_dir}/${sample_name}_stderr.log"

# Define the array of score_min parameters to test
score_min_params=(
    "L,0,-0.4"
    "L,0,-0.6"
    "L,0,-0.8"
    "L,0,-1.0"
    "L,-1,-0.6"
)

# Loop through each score_min parameter
for score_min in "${score_min_params[@]}"; do
    echo "Running Bismark for sample ${sample_name} with score_min ${score_min}"
    
    # Create a subdirectory for this parameter
    param_output_dir="${output_dir}/${sample_name}_score_${score_min//,/}"
    mkdir -p ${param_output_dir}

    # Run Bismark alignment
    bismark \
        -genome ${genome_folder} \
        -p 8 \
        -u 10000 \
        -score_min ${score_min} \
        -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
        -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
        -o ${param_output_dir} \
        --basename ${sample_name}_${score_min//,/} \
        2> "${param_output_dir}/${sample_name}-bismark_summary.txt"

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Sample ${sample_name} with score_min ${score_min} processed successfully."
    else
        echo "Sample ${sample_name} with score_min ${score_min} failed. Check ${stderr_log} for details."
    fi
done

# Mark the sample as completed in the checkpoint file
if [ $? -eq 0 ]; then
    echo ${sample_name} >> ${checkpoint_file}
    echo "All tests for sample ${sample_name} completed."
else
    echo "Sample ${sample_name} encountered errors. Check logs for details."
fi

# Define directories
summary_file="${output_dir}/parameter_comparison_summary.csv"

# Initialize summary file
echo "Sample,Score_Min,Alignment_Rate" > ${summary_file}

# Loop through parameter output directories
for dir in ${output_dir}/*_score_*; do
    if [ -d "$dir" ]; then
        # Extract sample name and score_min parameter from directory name
        sample_name=$(basename "$dir" | cut -d'_' -f1-4)
        score_min=$(basename "$dir" | grep -o "score_.*" | sed 's/score_//; s/_/,/g')

        # Locate the summary file
        summary_file_path="${dir}/${sample_name}_${score_min}_PE_report.txt"

        # Extract metrics
        mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{print "mapping efficiency ", $3}')
        
        # Append to the summary file
        echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
    fi
done
```

#### final script:

The best minscore for alignment for all samples was L,0,-1.0 

<img width="800" alt="minscoregraph" src="08-Bismark-Alignment-Assesment-images/alignment_bismark_minscore.png?raw=true">

However, the low mapping rates of higher quality samples (based on FastQC) is suspicious. Maybe those higher quality sequences are too long and need to be trimmed further to increase mapping?

```
nano scripts/bismark_align.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=100GB
#SBATCH -t 7-24:00:00
#SBATCH -q long #job lasting over 2 days
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

mkdir -p output_WGBS/align #make output directory if it does not exist

output_dir="output_WGBS/align"
checkpoint_file="output_WGBS/align/completed_samples.log"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}trimmed_LCM_*R1_001.fastq.gz)
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
    -genome ${genome_folder} \
    -p 8 \
    -score_min L,0,-1.0 \
    -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
    -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
    -o ${output_dir} \
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
    sample_name=$(basename "$file" | cut -d'_' -f1-4)
    score_min="L0-1.0"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
done
```

### Hmm. Interp of Bismark

Alignment rates are low, and are the lowest for libraries that had the highest initial quality. Could this be because shorter reads are aligning better?

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_raw_r1.png?raw=true" height="400">
<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_raw_r2.png?raw=true" height="400">
<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_trimmed_r1.png?raw=true" height="400">
<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_trimmed_r2.png?raw=true" height="400">

#### how to interpret this:

What does score-min mean?

https://github.com/FelixKrueger/Bismark/blob/f88314914725242c59289a534271b66cf9d461e3/docs/options/alignment.md: 

- "Sets a function governing the minimum alignment score needed for an alignment to be considered "valid" (i.e. good enough to report)."
- "Specifying L,0,-0.2 sets the minimum-score function f to f(x) = 0 + -0.2 * x, where x is the read length."

So, let's take score_min="L,0,-0.2" and "L,0,-1.0". What are different minimum score values for different read lengths?

- f(x) = 0 + -0.2 * (120) = -24
- f(x) = 0 + -0.2 * (130) = -26
- f(x) = 0 + -0.2 * (140) = -28

- f(x) = 0 + -1.0 * (120) = -120
- f(x) = 0 + -1.0 * (130) = -130
- f(x) = 0 + -1.0 * (140) = -140

So the longer the read, on average, the lower the minimum score. And by increasing the multiplier of X, we are decreasing the threshold for alignments to be considered valid. But, why would this decrease the mapping rate of samples with longer reads? Longer reads should be held to a lower cutoff given this formula.

## Plan: use --next-seq trimming to account for NovaSeq sequencing parameters:

https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types

"Quality trimming of reads using two-color chemistry (NextSeq)
Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) **encodes a G**. However, dark cycles also occur when sequencing “falls off” the end of the fragment. **The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.**

Since the regular quality-trimming algorithm cannot deal with this situation, you need to use the --nextseq-trim option:

**cutadapt --nextseq-trim=20 -o out.fastq input.fastq**
This works like regular quality trimming (where one would use -q 20 instead), except that the qualities of G bases are ignored."


## Trimming V2: Next-Seq trim

From the Zymo trio protocol:

"Libraries should be trimmed to remove any adapter sequence. No other trimming is required. Use the following sequences to trim the adapters:"

Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

I am going to use [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and quality control

Example code with comments:

```
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # Zymo Adapter Read1 
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # Zymo Adapter Read2
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \ #output files
    input_R1.fastq.gz input_R2.fastq.gz \ #input files
    -q 20,20 \ #trims low-quality bases (score < 20) from the 3' end (first 20) and 5' (second 20) of the read
    --nextseq-trim=20 \ #ignores the quality of G bases for trimming at the end of reads to account for NextSeq methods
    --minimum-length 20 #after trimming, only keep a sequence if longer than 20 bp
```

```
nano scripts/wgbs_cutadapt_V2.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

# load modules needed

module load uri/main
module load cutadapt/3.5-GCCcore-11.2.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' LCM*R1*.fastq.gz))
R2_raw=($('ls' LCM*R2*.fastq.gz))

for i in ${!R1_raw[@]}; do
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o trimmed_V2_${R1_raw[$i]} -p trimmed_V2_${R2_raw[$i]} \
    ${R1_raw[$i]} ${R2_raw[$i]} \
    -q 20,20 --nextseq-trim=20 --minimum-length 20 --cores=20

    echo "trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# unload conflicting modules with modules needed below
module unload cutadapt/3.5-GCCcore-11.2.0

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

#make trimmed_V2_qc output folder
mkdir ../output_WGBS/trimmed_V2_qc/

# Create an array of fastq files to process
files=($('ls' trimmed_V2*.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/trimmed_V2_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ../output_WGBS/trimmed_V2_qc/

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)

```

Okay, this made tiny improvements to the adapter content in the R2 reads, but did not have a major impact on QC otherwise.

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_WGBS/trimmed_V2_qc), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/72ae2716b6de1e90b76d4e522721a445dd47c04c/output_WGBS/trimmed_V2_qc/multiqc_report.html)

## Trimming V3: Next-Seq and increase minimum post-trim length

Remove reads that are less than 40 bp long. 

We have a little peak of reads that were < 40 bp post-trimming, and I want to see if removing these helps with alignment at all.

<img src="../output_WGBS/trimmed_qc/multiqc_screenshots/fastqc_sequence_length_distribution_plot.png?raw=true" height="400">

From the Zymo trio protocol:

"Libraries should be trimmed to remove any adapter sequence. No other trimming is required. Use the following sequences to trim the adapters:"

Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

I am going to use [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) for trimming and quality control

Example code with comments:

```
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # Zymo Adapter Read1 
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \  # Zymo Adapter Read2
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \ #output files
    input_R1.fastq.gz input_R2.fastq.gz \ #input files
    -q 20,20 \ #trims low-quality bases (score < 20) from the 3' end (first 20) and 5' (second 20) of the read
    --nextseq-trim=20 \ #ignores the quality of G bases for trimming at the end of reads to account for NextSeq methods
    --minimum-length 40 #after trimming, only keep a sequence if longer than 20 bp
```

```
nano scripts/wgbs_cutadapt_V3.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
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
module load cutadapt/3.5-GCCcore-11.2.0

#make arrays of R1 and R2 reads
R1_raw=($('ls' LCM*R1*.fastq.gz))
R2_raw=($('ls' LCM*R2*.fastq.gz))

for i in ${!R1_raw[@]}; do
    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o trimmed_V3_${R1_raw[$i]} -p trimmed_V3_${R2_raw[$i]} \
    ${R1_raw[$i]} ${R2_raw[$i]} \
    -q 20,20 --nextseq-trim=20 --minimum-length 40 --cores=20

    echo "trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# unload conflicting modules with modules needed below
module unload cutadapt/3.5-GCCcore-11.2.0

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load all/MultiQC/1.12-foss-2021b

#make trimmed_V3_qc output folder
mkdir ../output_WGBS/trimmed_V3_qc/

# Create an array of fastq files to process
files=($('ls' trimmed_V3*.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/trimmed_V3_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ../output_WGBS/trimmed_V3_qc/

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
```

### Trimming version 3 interpretation

- Trimming V2 (adding Next-Seq parameter) made basically no difference except slightly decreased adapter content. 
  - Interestingly, for both V2 and V3 when I added this parameter there is something funky going on at the end of the reads where the G content drops super low and it causes the per base sequence content to fail. 
  - Original Trimming:
    - <img src="../output_WGBS/trimmed_qc/multiqc_screenshots/PerBaseSeqContent_NextSeqEX.png?raw=true" height="400">
  - NextSeq Trimming:
    - <img src="../output_WGBS/trimmed_V3_qc/multiqc_screenshots/PerBaseSeqContent_NextSeqEX.png?raw=true" height="400">
  - **We'll go with it for now, because the V3 trimming worked so well. If there is an issue, I will keep going without the nextSeq trim and just increase the min_length filter to 40 as in V3 trim**
- Trimming V3, with the Next-Seq trimming and increased min-length filter worked really well. Basically no adapter content left. Going to align these now.

[View results here](https://github.com/zdellaert/LaserCoral/tree/main/output_WGBS/trimmed_V3_qc), [MultiQC report](https://github.com/zdellaert/LaserCoral/blob/72ae2716b6de1e90b76d4e522721a445dd47c04c/output_WGBS/trimmed_V3_qc/multiqc_report.html)

<img src="../output_WGBS/trimmed_V3_qc/multiqc_screenshots/fastqc_adapter_content_plot.png?raw=true" height="400">
<img src="../output_WGBS/trimmed_V3_qc/multiqc_screenshots/fastqc_sequence_length_distribution_plot.png?raw=true" height="400">
<img src="../output_WGBS/trimmed_V3_qc/multiqc_screenshots/fastqc-status-check-heatmap.png?raw=true" height="400">
<img src="../output_WGBS/trimmed_V3_qc/multiqc_screenshots/overrep_seq.png?raw=true">


## Bismark Alignment of V3 Trimmed Reads

```
nano scripts/bismark_align_V3.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=100GB
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

mkdir -p output_WGBS/align_V3 #make output directory if it does not exist

output_dir="output_WGBS/align_V3"
checkpoint_file="output_WGBS/align_V3/completed_samples.log"

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
    -genome ${genome_folder} \
    -p 8 \
    -score_min L,0,-1.0 \
    -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
    -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
    -o ${output_dir} \
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
    sample_name=$(basename "$file" | cut -d'_' -f3-4)
    score_min="L0-1.0"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
done
```

### Interp of Bismark

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_trimmed_v3_r1.png?raw=true" height="400">
<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_vs_qc_trimmed_v3_r2.png?raw=true" height="400">

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_trimmed_v3_batch1.png?raw=true" height="400">
<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_trimmed_v3_batch2.png?raw=true" height="400">

Okay , so an improvement but nothing super substantial. I am considering trimming all reads to length 100bp to try to get rid of any potential length biasing in the alignment? I could also try bwa-meth or hisat2 as alternative aligners.

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_orig_vs_v3.png?raw=true" height="400">

But, what I did not realize was that the next step in the Bismark pipeline, deduplication, can increase the effective alignment rate. Since I have high rates of duplication, I will try this before looking into alternative mapping or more aggressive trimming strategies.

## Deduplication

```
nano scripts/bismark_dedup.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load parallel/20240822
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

mkdir -p output_WGBS/dedup_V3 #make output directory if it does not exist

find output_WGBS/align_V3/*.bam | \
xargs -n 1 basename -s .bam | \
parallel -j 8 deduplicate_bismark \
--bam \
--paired \
--output_dir output_WGBS/dedup_V3 \
output_WGBS/align_V3/{}.bam
```

Oh no! About 87-95% of alignments were removed for each sample in this step, which decreased the alignment rate even furhter!

Let's just do the methylation extraction to see:

### Extract

```
nano scripts/bismark_extract.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b

cd output_WGBS/dedup_V3

bismark_methylation_extractor \
--bedGraph \
--counts \
--comprehensive \
--merge_non_CpG \
--multicore 28 \
--buffer_size 75% \
*deduplicated.bam
```

## Methylation call

```
nano scripts/bismark_call.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=50GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load parallel/20240822
module load uri/main
module load Bismark/0.23.1-foss-2021b

cd output_WGBS/dedup_V3

find *deduplicated.bismark.cov.gz | \
xargs -n 1 basename -s _pe.deduplicated.bismark.cov.gz | \
parallel -j 10 coverage2cytosine \
--genome_folder ../../references/ \
-o {} \
--merge_CpG \
--zero_based \
{}_pe.deduplicated.bismark.cov.gz
```

### Reports

Turns out everything has to be in the same folder:

```
salloc
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load all/MultiQC/1.12-foss-2021b

cd output_WGBS
mv output_WGBS/align_V3/* output_WGBS/dedup_V3/

cd output_WGBS/dedup_V3
bismark2report
bismark2summary *pe.bam

bam2nuc --genome_folder ../../references/ *_pe.deduplicated.bam

multiqc .
```



<img src="../output_WGBS/dedup_V3/report_screenshots/Summary_1.png?raw=true" height="400">
<img src="../output_WGBS/dedup_V3/report_screenshots/Summary_2.png?raw=true" height="400">

## Thoughts and next steps. 

Okay, the alignment is not great. And I don't like how the best looking libraries (32, 33) are aligning the worst. I am going to hard trim all reads to 100 bp and try to align these and see if this reduces possible length bias.

**Something weird is happening with the GC content.**

- I think i need to trim polyA. 


## I am running out of space. 

### Make a scratch directory

https://docs.unity.uri.edu/documentation/managing-files/hpc-workspace/

```
ws_allocate -G pi_hputnam_uri_edu -m zdellaert@uri.edu -r 2 shared 30
```

Info: creating workspace.
/scratch3/workspace/zdellaert_uri_edu-shared
remaining extensions  : 3
remaining time in days: 30


## Hard trim reads to 100 bp and see if this reduces length dependence

```
nano scripts/wgbs_trim_V4.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
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
to_trim=($('ls' trimmed_V3_LCM*.fastq.gz))

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

# go to directory with trimmed files

cd /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V4

# Create an array of fastq files to process
files=($('ls' *100bp_5prime.fq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ../output_WGBS/trimmed_V4_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ../output_WGBS/trimmed_V4_qc/

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
```

## Bismark Alignment of V4 Trimmed Reads

```
nano scripts/bismark_align_V4.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-9 #for 10 samples
#SBATCH --mem=100GB
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
reads_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V4/" #directory containing trimmed data to align
genome_folder="references/" #directory containing original unmodified genome fasta and bismark Bisulfite_Genome directory

mkdir -p output_WGBS/align_V4 #make output directory if it does not exist

output_dir="output_WGBS/align_V4"
checkpoint_file="output_WGBS/align_V4/completed_samples.log"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}*100bp_5prime.fq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.100bp_5prime.fq.gz")

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
    -genome ${genome_folder} \
    -p 8 \
    -score_min L,0,-1.0 \
    -1 ${reads_dir}${sample_name}_R1_001.100bp_5prime.fq.gz \
    -2 ${reads_dir}${sample_name}_R2_001.100bp_5prime.fq.gz \
    -o ${output_dir} \
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
    score_min="L0-1.0"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${short_name},${score_min},${mapping}" >> ${summary_file}
done
```

## Interpretation

This did not increase alignment rates, but roughly decreased the alignment rate of each sample by ~1%.

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_orig_vs_v4.png?raw=true" height="400">

---


## Side project: get methylseq to work:

https://nf-co.re/methylseq/3.0.0/

"ANNOUNCEMENTS 📢
💡 Unity profile is now available for Nextflow nf-core pipelines
An institutional Unity profile is now available for Nextflow nf-core pipelines (docs, config)!

To use it, add the -profile unity option to any nf-core pipeline. This allows for each process to be submitted as a slurm job and use apptainer for dependency management."

```
nano samplesheet.csv
```

paste this:

sample,fastq_1,fastq_2,genome
SRR389222_sub1,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz,,
SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz,,
SRR389222_sub3,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz,,

```
nano scripts/methylseq_test.sh 
```

test script:

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

## Load Nextflow and Apptainer environment modules
module purge
module load nextflow/24.10.3
module load apptainer/latest

nextflow run nf-core/methylseq \
  --input samplesheet.csv \
  --genome GRCh38 \
  -profile test,unity
```

### Run with default parameters

```
nano data_WGBS/LCM_methylseq_input.csv
```

```
sample,fastq_1,fastq_2,genome
LCM_1,data_WGBS/LCM_1_S1_R1_001.fastq.gz,data_WGBS/LCM_1_S1_R2_001.fastq.gz,
LCM_3,data_WGBS/LCM_3_S2_R1_001.fastq.gz,data_WGBS/LCM_3_S2_R2_001.fastq.gz,
LCM_11,data_WGBS/LCM_11_S3_R1_001.fastq.gz,data_WGBS/LCM_11_S3_R2_001.fastq.gz,
LCM_12,data_WGBS/LCM_12_S4_R1_001.fastq.gz,data_WGBS/LCM_12_S4_R2_001.fastq.gz,
LCM_17,data_WGBS/LCM_17_S7_R1_001.fastq.gz,data_WGBS/LCM_17_S7_R2_001.fastq.gz,
LCM_18,data_WGBS/LCM_18_S8_R1_001.fastq.gz,data_WGBS/LCM_18_S8_R2_001.fastq.gz,
LCM_24,data_WGBS/LCM_24_S5_R1_001.fastq.gz,data_WGBS/LCM_24_S5_R2_001.fastq.gz,
LCM_25,data_WGBS/LCM_25_S6_R1_001.fastq.gz,data_WGBS/LCM_25_S6_R2_001.fastq.gz,
LCM_32,data_WGBS/LCM_32_S9_R1_001.fastq.gz,data_WGBS/LCM_32_S9_R2_001.fastq.gz,
LCM_33,data_WGBS/LCM_33_S10_R1_001.fastq.gz,data_WGBS/LCM_33_S10_R2_001.fastq.gz,
```

```
nano scripts/methylseq.sh 
```


```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH -t 48:00:00
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

#nextflow run nf-core/methylseq \
#  --input ./data_WGBS/LCM_methylseq_input.csv \
#  --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_default_test \
#  --igenomes_ignore \
#  --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
#  --relax_mismatches \
#  -profile unity \
#  -name methylseq_default_test

nextflow run nf-core/methylseq \
    --input ./data_WGBS/LCM_methylseq_input.csv \
    --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_default_test \
    --igenomes_ignore \
    --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
    --relax_mismatches \
    -profile unity \
    -resume
```


Emma additional params:

--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--non_directional \
--cytosine_report \
--unmapped \

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

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_bowtie_vs_hisat.png?raw=true" height="400">


## Trimming V5: Poly A

```
nano scripts/wgbs_trim_V5.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
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
```

### Trimming version 5 interpretation

This got rid of poly-A but the GC content graphs are still double-peaked enough to the point where I'm not going to try aligning these yet at this stage.

## Trimming V6: FastP

Inspired by Sam White's work here:
https://robertslab.github.io/sams-notebook/posts/2025/2025-01-02-Trimming---A.pulchra-WGBS-with-fastp-FastQC-and-MultiQC-on-Raven/

```
nano scripts/wgbs_trim_V6.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
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
  2> ../scripts/outs_errs/"$sample_name".fastp-trim.stderr

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
files=($('ls' trimmed_V6*.fastq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

cd ${qc_dir}

#Compile MultiQC report from FastQC files
multiqc *  #Compile MultiQC report from FastQC files 

echo "QC of trimmed data complete." $(date)
```

This did not work well at all. 

## Trimming V7: Flexbar

Install:

```
cd work/pi_hputnam_uri_edu/pgrams

wget https://github.com/seqan/flexbar/releases/download/v3.5.0/flexbar-3.5.0-linux.tar.gz
tar xzf flexbar-3.5.0-linux.tar.gz

echo 'export PATH=/work/pi_hputnam_uri_edu/pgrams/flexbar-3.5.0-linux:$PATH' >> ~/.bash_profile
echo 'export LD_LIBRARY_PATH=/work/pi_hputnam_uri_edu/pgrams/flexbar-3.5.0-linux:$LD_LIBRARY_PATH' >> ~/.bash_profile

source ~/.bash_profile

flexbar -h
```

```
nano scripts/wgbs_trim_V7_test2.sh
```

```
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

base_name=${R1_raw[$i]%%_R1*.fastq.gz}

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
    --threads 120 \
    -t ${base_name}_fbtrim

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

for i in ${!R1_raw[@]}; do

base_name=${R1_raw[$i]%%_R1*.fastq.gz}

    flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    -a adapters.fasta \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \
    --zip-output GZ \
    --threads 120 \
    -t ${base_name}_fbtrim_test_polyAdapt

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
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-left AG --htrim-right AG \
    --zip-output GZ \
    --threads 120 \
    -t ${base_name}_fbtrim_test_polyAdaptTrim

    echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

for i in ${!R1_raw[@]}; do

base_name=${R1_raw[$i]%%_R1*.fastq.gz}

    flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim WIN --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \
    --zip-output GZ \
    --threads 120 \
    -t ${base_name}_fbtrim_test_WIN

    echo "Trimming of ${R1_raw[$i]} and ${R2_raw[$i]} complete"
done

# load modules needed
module purge
module load parallel/20240822
module load fastqc/0.12.1

# Create an array of fastq files to process
files=($('ls' *fbtrim_test*gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)
```

Example commented code:

```
    flexbar \
    -r ${R1_raw[$i]} \
    -p ${R2_raw[$i]} \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \ #could decrease this
    --adapter-error-rate 0.1 \ #could increase this
    --adapter-trim-end RIGHT \ #could change to ANY
    --adapter-pair-overlap ON \ #default is off, but recommended to be on
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right G \ #trim poly-G
    --zip-output GZ \
    --threads 20 \
    -t ${base_name}_fbtrim
```

### assessment of flexbar.

I ran a lot of tests on sample LCM_1 and assessed the results with FastQC. 

Changes from default are listed in itialcs. Additional parameters added are in bold.

1. "Default"
    - max-uncalled:          0
    - min-read-length:       40
    - qtrim:                 TAIL
    - qtrim-format:          sanger
    - qtrim-threshold:       20  (53)
    - htrim-right:           G
    - htrim-min-length:      3
    - htrim-error-rate:      0.1
    - adapter-pair-overlap:  ON
    - adapter-trim-end:      RIGHT
    - adapter-min-overlap:   3
    - adapter-min-poverlap:  40
    - adapter-error-rate:    0.1
    - adapter-match:         1
    - adapter-mismatch:     -1
    - adapter-gap:          -6
    - Adapter: TruSeq preset
      - AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    - Adapter2: TruSeq preset
      - AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
4. "WIN" = Sliding window quality trimming
   - *qtrim:                 WIN*
   - **qtrim-win-size:        5**
2. "polyAdapt": Add polyA and polyG as "adapter sequences"
   - *Adapters provided in fasta file:*
     - adapter1               AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
     - adapter2               AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
     - polyA                  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     - polyG                  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
3. "poltAdaptTrim": Add polyA and polyG as adapters, also do additional polyA and polyG trimming
   - *Adapters provided in fasta file:*
     - adapter1               AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
     - adapter2               AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
     - polyA                  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     - polyG                  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
   - **htrim-left:            AG**
   - *htrim-right:           AG*
4. "ANY": Look for internal adapters 
   - *Adapters provided in fasta file:*
     - adapter1               AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
     - adapter2               AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
     - polyA                  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     - polyG                  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
   - *adapter-trim-end:      ANY*
   - ~~*adapter-pair-overlap:  ON*~~
   - ~~*adapter-min-poverlap:  40*~~
5. "RLG": Internal adapters, plus lefthand trimming of polyG
   - *Adapters provided in fasta file:*
     - adapter1               AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
     - adapter2               AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
     - polyA                  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
     - polyG                  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
   - **htrim-left:           G**
   - *adapter-trim-end:      ANY*
   - ~~*adapter-pair-overlap:  ON*~~
   - ~~*adapter-min-poverlap:  40*~~

Quick assessment: 
1. default
   1. worked very similarly to cutadapt
   2. Most polyA and polyG of all tests
2. sliding window
   1. improved overall quality, decreased final number of seqeunces kept
   2. decreased adapter content compared to default
   3. still a lot of polyA
3. PolyAdapt
   1. Higher adapter content than default
4. PolyAdapt with trimming
   1. Higher adapter content than default
   2. Had least polyA and polyG at end
5. ANY
   1. Higher adapter content than default
   2. did really odd things to the reads
   3. not recommended
6. RLG
   1. Higher adapter content than default
   2. did really odd things to the reads
   3. not recommended

Based on the QC of these tests, I am going to employ the following changes to the script:

- Sliding window trimming resulted in the best overall quality of reads and per base sequent content
  - **Perform window quality trimming, consider testing window size parameter**
- For all reads, the first 10 bp is still shaky in per base sequent content
  - Consider hard trim of first 10 bp for all samples and reads
- Consider decreasing final read length cutoff to maintain more reads
- Adding polyA and polyG as adapters decreased removal of the truseq adapters
- *Is trimming polyA appropriate if we do not know why it is there?*
- **GC content issue is still not solved.**
- 

1. Sliding window trimming
   1. *qtrim:                 WIN*
   2. **qtrim-win-size:        5**
2. Add right and lefthand polyG trimming
   1. -htrim-right: G -htrim-left: G


## Figuring out GC content issues

1. let's do bam 2 qc on some of the bismark alignments and then look at the GC content of the reads that successfully aligned

```
nano scripts/bam2fq.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=40
#SBATCH --mem=500GB
#SBATCH -t 48:00:00
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
```

Okay. I looked at these in geneious and there are Many sequences with >90% C, T, or A. I am going to try a series of polyX trims using LCM_1 as my test subject:

## Flexbar Test "3"

```
nano scripts/wgbs_trim_V7_test3.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=16 
#SBATCH --mem=100GB 
#SBATCH -t 06:00:00 
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%A_%a"
#SBATCH --output=../scripts/outs_errs/"%x_output.%A_%a"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
#SBATCH --array=0-6

module load fastqc/0.12.1

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/further_flexbar"
mkdir -p $output_dir

# Define trimming variations
trim_types=("trim_G" "trim_5G" "trim_T" "trim_C" "trim_A" "trim_all_3g" "trim_all")

trim_cmds=(
    "--htrim-right G --htrim-left G"
    "--htrim-right G --htrim-left G --htrim-min-length 5"
    "--htrim-right T --htrim-left T --htrim-min-length 5"
    "--htrim-right C --htrim-left C --htrim-min-length 5"
    "--htrim-right A --htrim-left A --htrim-min-length 5"
    "--htrim-right GCTA --htrim-left GCTA --htrim-min-length2 5"
    "--htrim-right GCTA --htrim-left GCTA --htrim-min-length 5"
)

# Get the trimming command for this array task
trim_type=${trim_types[$SLURM_ARRAY_TASK_ID]}
trim_cmd=${trim_cmds[$SLURM_ARRAY_TASK_ID]}

echo "Running Flexbar with ${trim_type} trimming..."

flexbar \
    -r LCM_1_S1_R1_001.fastq.gz \
    -p LCM_1_S1_R2_001.fastq.gz \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    ${trim_cmd} \
    --zip-output GZ \
    --threads 16 \
    -t "${output_dir}/${trim_type}"

# Run FastQC on the output
fastqc "${output_dir}/${trim_type}*.fastq.gz" -o $output_dir

echo "Finished processing ${trim_type}"
```

### Based on QC, the "Trim All" worked the best, all with a min length of 5:

```
flexbar \
    -r LCM_1_S1_R1_001.fastq.gz \
    -p LCM_1_S1_R2_001.fastq.gz \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right GCTA --htrim-left GCTA --htrim-min-length 5 \
    --zip-output GZ \
    --threads 16 \
    -t "${output_dir}/${trim_type}"
```

Buuuuut the GC content is barely changed.

Further things to try: 
- -u, --max-uncalled  Allowed uncalled bases N for each read. Default: 0.
  - **could increase this.**
- decrease min read length back to 20ish
- -he, --htrim-error-rate  Error rate threshold for mismatches. Default: 0.1.
  - **could increase this**

```
nano scripts/wgbs_trim_V7_test4.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=16 
#SBATCH --mem=100GB 
#SBATCH -t 06:00:00 
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%A_%a"
#SBATCH --output=../scripts/outs_errs/"%x_output.%A_%a"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
#SBATCH --array=0-5

module load fastqc/0.12.1

output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/further_flexbar"
mkdir -p $output_dir

# Define argument variations
arg_types=("trim_all_plus" "trim_all_plus_lenient" "trim_all_plus_more_lenient" "trim_all_plus_lenient_max10" "trim_all_plus_lenient_qpost" "trim_all_plus_lenient_pretrim")

cmds=(
    "--max-uncalled 2"
    "--max-uncalled 2 --htrim-error-rate 0.2"
    "--max-uncalled 2 --htrim-error-rate 0.4"
    "--max-uncalled 2 --htrim-error-rate 0.2 -htrim-min-length2 10"
    "--max-uncalled 2 --htrim-error-rate 0.2 --qtrim-post-removal"
    "--max-uncalled 2 --htrim-error-rate 0.2 --pre-trim-left 10"
)

# Get the trimming command for this array task
arg_type=${arg_types[$SLURM_ARRAY_TASK_ID]}
cmd=${cmds[$SLURM_ARRAY_TASK_ID]}

flexbar \
    -r LCM_1_S1_R1_001.fastq.gz \
    -p LCM_1_S1_R2_001.fastq.gz \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right GCTA --htrim-left GCTA --htrim-min-length 5 \
    ${cmd} \
    --zip-output GZ \
    --threads 16 \
    -t "${output_dir}/${arg_type}"

# Run FastQC on the output
fastqc "${output_dir}/${arg_type}_1.fastq.gz" -o $output_dir
fastqc "${output_dir}/${arg_type}_2.fastq.gz" -o $output_dir

echo "Finished processing ${arg_type}"
```

### this worked the best:

flexbar \
    -r LCM_1_S1_R1_001.fastq.gz \
    -p LCM_1_S1_R2_001.fastq.gz \
    --adapter-preset TruSeq \
    --adapter-min-overlap 3 \
    --adapter-error-rate 0.1 \
    --adapter-trim-end RIGHT \
    --adapter-pair-overlap ON \
    --qtrim TAIL --qtrim-threshold 20 \
    --qtrim-format sanger \
    --min-read-length 40 \
    --htrim-right GCTA --htrim-left GCTA --htrim-min-length 5 \
    --max-uncalled 2 --htrim-error-rate 0.2 -htrim-min-length2 10 \
    --zip-output GZ \
    --threads 16 \
    -t "${output_dir}/${arg_type}"

## Okay -- None of these polyX trims are getting rid of the huge GC content peak.

Let's try a low complexity filter.



## Trimming V8: FastP with low complexity filter

Initial FastP script (V6) was inspired by Sam White's work here:
https://robertslab.github.io/sams-notebook/posts/2025/2025-01-02-Trimming---A.pulchra-WGBS-with-fastp-FastQC-and-MultiQC-on-Raven/

```
nano scripts/wgbs_trim_V8_test.sh #write script for first trimming pass and QC, enter text in next code chunk
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=16 
#SBATCH --mem=100GB 
#SBATCH -t 06:00:00 
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%A_%a"
#SBATCH --output=../scripts/outs_errs/"%x_output.%A_%a"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS
#SBATCH --array=0-5

# load modules needed

module load uri/main
module load fastp/0.23.2-GCC-11.2.0

#make trimmed_V8_qc output folder
mkdir -p /home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V8_qc
mkdir -p /scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V8/
output_dir="/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V8/"
qc_dir="/home/zdellaert_uri_edu/LaserCoral/output_WGBS/trimmed_V8_qc"

# Define argument variations
arg_types=("default" "LowC_30" "LowC_15" "LowC_50" "LowC_70" "LowC_80")

cmds=(
    "--length_required 40"
    "--length_required 40 --low_complexity_filter"
    "--length_required 40 --low_complexity_filter --complexity_threshold 15"
    "--length_required 40 --low_complexity_filter --complexity_threshold 50"
    "--length_required 40 --low_complexity_filter --complexity_threshold 70"
    "--length_required 40 --low_complexity_filter --complexity_threshold 80"
)

# Get the trimming command for this array task
arg_type=${arg_types[$SLURM_ARRAY_TASK_ID]}
cmd=${cmds[$SLURM_ARRAY_TASK_ID]}

fastp \
  --in1 LCM_1_S1_R1_001.fastq.gz \
  --in2 LCM_1_S1_R2_001.fastq.gz \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --poly_g_min_len 3 \
  --trim_poly_x \
  --cut_front --cut_front_mean_quality 20 \
  ${cmd} \
  --thread 16 \
  --html ${qc_dir}/LCM_1_"${arg_type}"_fastp.report.html \
  --json ${qc_dir}/LCM_1_"${arg_type}"_fastp.report.json \
  --out1 ${output_dir}/LCM_1_R1_"${arg_type}"_fastp.fq.gz \
  --out2 ${output_dir}/LCM_1_R2_"${arg_type}"_fastp.fq.gz

    echo "trimming of LCM_1 with "${arg_type}" complete"

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Create an array of fastq files to process
files=($('ls' ${output_dir}/*_"${arg_type}"_fastp.fq.gz)) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ${qc_dir} && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

echo "QC of "${arg_type}"-trimmed data complete." $(date)
```



## Trimming V9?: BBtools

Install:

```
cd work/pi_hputnam_uri_edu/pgrams

wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz

tar -xvzf BBMap_39.01.tar.gz
cd bbmap
chmod +x *.sh

echo 'export PATH=/work/pi_hputnam_uri_edu/pgrams/bbmap:$PATH' >> ~/.bash_profile
source ~/.bash_profile

bbduk.sh --help
```



```
nano scripts/wgbs_trim_V9_test.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=12
#SBATCH --mem=50GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/data_WGBS

source ~/.bash_profile

bbduk.sh in=trimmed_V3_LCM_1_S1_R1_001.fastq.gz out=trimmed_V3_LCM_1_S1_R1_001_filtered.fastq.gz outm=trimmed_V3_LCM_1_S1_R1_001_removed.fastq.gz entropy=0.5

R1_IN="trim_all_plus_lenient_max10_1.fastq.gz"
R2_IN="trim_all_plus_lenient_max10_2.fastq.gz"
R1_OUT="trim_all_plus_lenient_max10_1_bbfiltered.fastq.gz"
R2_OUT="trim_all_plus_lenient_max10_2_bbfiltered.fastq.gz"
R1_REMOVED="trim_all_plus_lenient_max10_1_bbremoved.fastq.gz"
R2_REMOVED="trim_all_plus_lenient_max10_2_bbremoved.fastq.gz"

cd /home/zdellaert_uri_edu/zdellaert_uri_edu-shared/data_WGBS/further_flexbar/

bbduk.sh in1=${R1_IN} in2=${R2_IN} out1=${R1_OUT} out2=${R2_OUT} outm1=${R1_REMOVED} outm2=${R2_REMOVED} entropy=0.5

R1_IN="trim_all_plus_lenient_max10_1.fastq.gz"
R2_IN="trim_all_plus_lenient_max10_2.fastq.gz"
R1_OUT="trim_all_plus_lenient_max10_1_bbfiltered_E75.fastq.gz"
R2_OUT="trim_all_plus_lenient_max10_2_bbfiltered_E75.fastq.gz"
R1_REMOVED="trim_all_plus_lenient_max10_1_bbremoved_E75.fastq.gz"
R2_REMOVED="trim_all_plus_lenient_max10_2_bbremoved_E75.fastq.gz"

bbduk.sh in1=${R1_IN} in2=${R2_IN} out1=${R1_OUT} out2=${R2_OUT} outm1=${R1_REMOVED} outm2=${R2_REMOVED} entropy=0.75

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1

# Create an array of fastq files to process
files=(ls *bb*) 

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 20 "fastqc {} -o ./ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

echo "QC of "${arg_type}"-trimmed data complete." $(date)
```