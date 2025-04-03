# LCM WGBS Bioinformatic Processing

Script Written By:  Dellaert
Last Updated: 2/24/2025

- Sample prep: https://zdellaert.github.io/ZD_Putnam_Lab_Notebook/LCM-Exp-Sample-Prep/
- DNA extractions: 
- WGBS library prep: 

## Sequencing information

- Sequenced through Genohub Service Provider: Oklahoma Medical Research Foundation NGS Core
- Instrument: Illumina NovaSeq X Plus - 25B - PE 150 Cycle
- Read length: 2 x 150bp (Paired End)
- Number of samples: 10
- Guaranteed number of pass filter PE reads/sample: 400M (200M in each direction)

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

## Bismark Alignment to P. acuta genome: Attempt 1

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


## Trimming V2: Next-Seq trim

To try to improve alignment results, I am going to try alternative trimming parameters. I will use --next-seq trimming to account for NovaSeq sequencing parameters:

https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types

"Quality trimming of reads using two-color chemistry (NextSeq)
Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) **encodes a G**. However, dark cycles also occur when sequencing “falls off” the end of the fragment. **The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.**

Since the regular quality-trimming algorithm cannot deal with this situation, you need to use the --nextseq-trim option:

**cutadapt --nextseq-trim=20 -o out.fastq input.fastq**
This works like regular quality trimming (where one would use -q 20 instead), except that the qualities of G bases are ignored."

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


## Bismark Alignment: Attempt 2

I will now use bismark to align the V3-trimmed reads to the *P. acuta* genome.

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
- note, I tried the 100 bp test and it decreased alignment rates. same for testing hisat as an aligner. see `code/notes` for more info and code used.

<img src="08-Bismark-Alignment-Assesment-images/alignment_bismark_orig_vs_v3.png?raw=true" height="400">

But, what I did not realize was that the next step in the Bismark pipeline, deduplication, can increase the effective alignment rate. Since I have high rates of duplication, I will try this before looking into alternative mapping or more aggressive trimming strategies.

### Deduplication

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

Oh no! About 87-95% of alignments were removed for each sample in this step, which decreased the alignment rate even further!

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

### Methylation call

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

### Sort deduplicated bams

```
nano scripts/bismark_sort_dedup.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load samtools/1.19.2

cd output_WGBS/dedup_V3

for bamFile in *deduplicated.bam; do
	prefix=$(basename $bamFile .bam)

	samtools cat ${bamFile} | samtools sort --threads 48 -o ${prefix}.sorted.bam
done
```

### Sort .cov files

from here on this is based on https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2021-10-21-KBay-Bleaching-Pairs-WGBS-Analysis-Pipeline.md#filter_cov

```
nano scripts/bismark_sort_cov.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

cd output_WGBS/dedup_V3

for file in *merged_CpG_evidence.cov
do
  sample=$(basename "${file}" .CpG_report.merged_CpG_evidence.cov)
  bedtools sort -i "${file}" \
  > "${sample}"_sorted.cov
done
```

### Create bedgraph files

```
nano scripts/bismark_bed.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# Bedgraphs for 5X coverage 

cd output_WGBS/dedup_V3

for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${sample}"_5x_sorted.bedgraph
done

# Bedgraphs for 10X coverage 

for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${sample}"_10x_sorted.bedgraph
done
```

### filter for 5x coverage

```
nano scripts/bismark_5xcov.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

cd output_WGBS/dedup_V3

### Filtering for CpG for 5x coverage. To change the coverage, replace X with your desired coverage in ($5+6 >= X)

for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${sample}"_5x_sorted.tab
done
```

```
# wc -l *5x_sorted.tab
   1297759 trimmed_V3_LCM_11_S3_5x_sorted.tab
   1563694 trimmed_V3_LCM_12_S4_5x_sorted.tab
   1982902 trimmed_V3_LCM_17_S7_5x_sorted.tab
   3275522 trimmed_V3_LCM_18_S8_5x_sorted.tab
   1850258 trimmed_V3_LCM_1_S1_5x_sorted.tab
   6722970 trimmed_V3_LCM_24_S5_5x_sorted.tab
   3813379 trimmed_V3_LCM_25_S6_5x_sorted.tab
   1323643 trimmed_V3_LCM_32_S9_5x_sorted.tab
    414843 trimmed_V3_LCM_33_S10_5x_sorted.tab
   1296247 trimmed_V3_LCM_3_S2_5x_sorted.tab
  23541217 total
```

### intersection of all samples

```
nano scripts/bismark_intersect.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

cd output_WGBS/dedup_V3

multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

#change number after == to your number of samples

cat CpG.all.samps.5x_sorted.bed | awk '$4 ==10' > CpG.filt.all.samps.5x_sorted.bed 
```

### intersection of all samples with transcripts (Gene bodies - GBM)

```
nano scripts/bismark_intersect_bed.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

awk '{if ($3 == "transcript") {print}}' references/Pocillopora_acuta_HIv2.gtf  > references/Pocillopora_acuta_HIv2_transcripts.gtf

cd output_WGBS/dedup_V3

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b ../../references/Pocillopora_acuta_HIv2_transcripts.gtf \
  > ${i}_gene
done
```

### keep only loci found in all samples

```
nano scripts/bismark_intersect_bed_all.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

awk '{if ($3 == "transcript") {print}}' references/Pocillopora_acuta_HIv2.gtf  > references/Pocillopora_acuta_HIv2_transcripts.gtf

cd output_WGBS/dedup_V3

for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```

### alternative interesection script, with bedtools map

```
nano scripts/bismark_intersect_bed2.sh 
```


```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --mem=250GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load bedtools2/2.31.1

cd output_WGBS/dedup_V3

for i in *5x_sorted.bedgraph
do
   sample=$(basename ${i} _5x_sorted.bedgraph)
   bedtools map -a ../../references/Pocillopora_acuta_HIv2.gtf.cleaned.bed -b  ${i} -c 4 -o mean >  ${sample}_gene_body_methylation.txt

done
```

### Run qualimap on sorted deduplicated bams

```
nano scripts/bismark_qualimap.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=1 #split one task over multiple CPU
#SBATCH --mem=32GB
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load qualimap/2.2.1

cd output_WGBS/dedup_V3

mkdir -p qualimap/bamqc

for bamFile in *deduplicated.sorted.bam; do
	prefix=$(basename $bamFile .bam)

	qualimap \
    --java-mem-size=29491M \
    bamqc \
     \
    -bam ${bamFile}  \
     \
    -p non-strand-specific \
    --collect-overlap-pairs \
    -outdir qualimap/bamqc/${prefix} \
    -nt 6
done
```


### Subset bams for viewing in IGV

```
nano scripts/bismark_subset.sh 
```

```
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=1 #split one task over multiple CPU
#SBATCH --mem=32GB
#SBATCH -t 02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load samtools/1.19.2

cd output_WGBS/dedup_V3

mkdir -p dedup_bam_subsets

for bamFile in *deduplicated.sorted.bam; do
        prefix=$(basename $bamFile .bam .bam | cut -d'_' -f1-4)

        samtools view -b $bamFile Pocillopora_acuta_HIv2___Sc0000000:1-1000000 > ${prefix}.deduplicated.sorted.sub.bam
        samtools index "${prefix}.deduplicated.sorted.sub.bam" "${prefix}.deduplicated.sorted.sub.bam.bai"
done

mv *sub.bam* subsets/
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

<img src="../output_WGBS/dedup_V3/report_screenshots/multiqc_1.png?raw=true"><img src="../output_WGBS/dedup_V3/report_screenshots/multiqc_2.png?raw=true" height="400">

<img src="../output_WGBS/dedup_V3/report_screenshots/qualimap_coverage_histogram.png?raw=true" height="400">
<img src="../output_WGBS/dedup_V3/report_screenshots/qualimap_coverage_histogram_2.png?raw=true" height="400">

<img src="../output_WGBS/dedup_V3/report_screenshots/qualimap_coverage_cumulative.png?raw=true" height="400">
<img src="../output_WGBS/dedup_V3/report_screenshots/qualimap_insert_size.png?raw=true" height="400">

<img src="../output_WGBS/dedup_V3/report_screenshots/qualimap_gc_content.png?raw=true" height="400">

<img src="../output_WGBS/dedup_V3/report_screenshots/bismark_alignment.png?raw=true" height="400">
<img src="../output_WGBS/dedup_V3/report_screenshots/bismark_deduplication.png?raw=true" height="400">

## Thoughts and next steps. 

Okay, the alignment is not great. And I don't like how the best looking libraries (32, 33) are aligning the worst. The aligner is maybe having trouble with repetitive regions? high-GC content sequences are aligning very poorly.


## BWAmeth Alignment using Methylseq

I tested a lot of tools to minimize adapter and polyG content. See [here](https://github.com/zdellaert/LaserCoral/tree/main/code/notes/Trim_Saga.md). In the end, none of them improved parameters enough beyond the cutadapt trimming to justify using extra tools and steps.

I also tested hisat-2 as an alternative aligner, which showed no improvements from bowtie2 (both used with bismark). Finally, I tried BWA-meth, wrapped by [methylseq](https://nf-co.re/methylseq/3.0.0/). This greatly improved coverage for all samples, even though there was still discrepancy across samples.

> "ANNOUNCEMENTS 📢
💡 Unity profile is now available for Nextflow nf-core pipelines
An institutional Unity profile is now available for Nextflow nf-core pipelines (docs, config)!

> To use it, add the -profile unity option to any nf-core pipeline. This allows for each process to be submitted as a slurm job and use apptainer for dependency management."



```
nano data_WGBS/LCM_methylseq_input_V3bwa.csv
```

```
sample,fastq_1,fastq_2,genome
LCM_1,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_1_S1_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_1_S1_R2_001.fastq.gz,
LCM_3,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_3_S2_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_3_S2_R2_001.fastq.gz,
LCM_11,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_11_S3_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_11_S3_R2_001.fastq.gz,
LCM_12,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_12_S4_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_12_S4_R2_001.fastq.gz,
LCM_17,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_17_S7_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_17_S7_R2_001.fastq.gz,
LCM_18,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_18_S8_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_18_S8_R2_001.fastq.gz,
LCM_24,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_24_S5_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_24_S5_R2_001.fastq.gz,
LCM_25,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_25_S6_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_25_S6_R2_001.fastq.gz,
LCM_32,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_32_S9_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_32_S9_R2_001.fastq.gz,
LCM_33,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_33_S10_R1_001.fastq.gz,/scratch3/workspace/zdellaert_uri_edu-shared/data_WGBS/trimmed_V3/trimmed_V3_LCM_33_S10_R2_001.fastq.gz,
```


```
nano scripts/methylseq_V3_bwa.sh 
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --mem=100GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_80 #email you when job starts, stops and/or fails
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
   --input ./data_WGBS/LCM_methylseq_input_V3bwa.csv \
   --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test \
   --aligner bwameth \
   --methyl_kit \
   --igenomes_ignore \
   --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
   --save_reference \
   --run_qualimap \
   --save_align_intermeds \
   --skip_trimming \
   -profile unity \
   -name methylseq_V3_bwa_test
```

### Additional parameters to consider

parameters I did not use which Putnam lab people have used in the past:

--clip_r1 10 \
--clip_r2 10 \
--three_prime_clip_r1 10 --three_prime_clip_r2 10 \
--cytosine_report \
--unmapped \

### Reflections

Wow, I had a lot of issues. It was a saga. Methylseq kept quitting right after the bwa alignments, which took ~24 hours. I kept trying different ways of resuming the pipeline, but it would always start the alignments over. *This is a known issue: see [here, v3.0.0 - [2024-12-16]](https://github.com/nf-core/methylseq/blob/master/CHANGELOG.md)*. I tried many ways to try to skip the alignment step, but all with no luck. I was able to find 2 solutions, which thankfully produced identical results:

1) Run all post-bwameth alignment steps **Manually**. This was a pain! But looked through methylseq source code and was able to write code to run samtools, picard, methyldackel, and multiqc with the same parameters that methylseq uses 
   1) NOTE: I did this before figuring out how to trick methylseq into skipping the alignment, and after confirming the results were the same, I removed the results files from this version of the analysis from this repository to avoid confusion. However, I am saving the script just in case (see below). 

2) Trick methylseq! After a run where it has successfully completed the alignments and then given up (why, methylseq, why?), you can trick methylseq into using the bam files it already made by changing the following line of code in wherever your nextflow directory is:
   1) For me, the full path is: **`/home/zdellaert_uri_edu/.nextflow/assets/nf-core/methylseq/modules/nf-core/bwameth/align/main.nf`**
   2) You need to replace from "script" to "stub" with the following code, making adjustments so the path points to the output path of your methylseq run you are trying to resume. For me, that was **`/scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/alignments/`**

  ```
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_file = "${prefix}.bam"
    
    """
    ln -sf "/scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/alignments/${prefix}.bam" .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam_source: "pre-existing"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
  ```


   2) Then, rerun your methylseq run with the `-resume` flag at the end of the code (i replace the `-name XXX` flag with `-resume`)
   3) NOTE: this is a known bug. If they fix the bug in future versions of methylseq, you should be able to resume the pipeline when it crashes without issue.
      1) *This is a known issue: see [here, v3.0.0 - [2024-12-16]](https://github.com/nf-core/methylseq/blob/master/CHANGELOG.md)*
      2) "Note: bwameth/align module still needs fixing for not resuming from cache. So, its cache has been made lenient (Minimal input file metadata (name and size only) are included in the cache keys) in its config. This strategy provides a workaround for caching invalidation by current bwameth/align module requirement to touch the index files before alignment. An issue we hope to have fixed in a release soon."


### Methylseq BWA-post processing script (if methylseq doesn't want to work)

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=48
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"
#SBATCH -D /home/zdellaert_uri_edu/zdellaert_uri_edu-shared/methylseq_V3_bwa_test_back/bwameth/alignments

## This is an alternative script to run the bwa-meth post processing manually if methylseq quits after the alignment and doesn't want to resume....

module load samtools/1.19.2
module load apptainer/latest

reference="Pocillopora_acuta_HIv2.assembly.fasta"

for bamFile in *.bam; do
	prefix=$(basename $bamFile .bam)

	samtools cat ${bamFile} | samtools sort --threads 48 ${reference} -o ${prefix}.sorted.bam
done

for bamFile in *sorted.bam; do
	prefix=$(basename $bamFile .bam)

	samtools index --threads 1 $bamFile
    samtools flagstat --threads 48 $bamFile > ${prefix}.flagstat
	samtools stats --threads 48    $bamFile > ${prefix}.stats

	singularity exec https://depot.galaxyproject.org/singularity/picard:3.3.0--hdfd78af_0 picard MarkDuplicates \
	    --ASSUME_SORTED true --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT \
		--INPUT $bamFile \
    	--OUTPUT $bamFile.dedup.bam \
		--REFERENCE_SEQUENCE $reference \
    	--METRICS_FILE $prefix.MarkDuplicates.metrics.txt

	samtools index --threads 1 $bamFile.dedup.bam

    #output methylkit files
	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --methylKit $reference $bamFile.dedup.bam

    #output bismark-style cytosine report
    singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --CHG --CHH --cytosine_report $reference $bamFile.dedup.bam

    #output default bedgraph format
	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 $reference $bamFile.dedup.bam

	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel mbias $reference $bamFile.dedup.bam mbias/$prefix --txt > ${prefix}.mbias.txt
done

echo "Collecting MultiQC input files..."

multiqc_input_files=(
	*.MarkDuplicates.metrics.txt
	*.flagstat
	*.stats
	*.bam.bai
	*.bedGraph
	*.mbias.txt
)

module load uri/main
module load all/MultiQC/1.12-foss-2021b


# Generate MultiQC report
echo "Generating MultiQC report..."
#multiqc ${multiqc_input_files[@]} -o multiqc_report
multiqc *
```

## "Final" (as of now) BWA-Meth aligned MultiQC Report: 

<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/multiqc_1.png?raw=true">

<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_coverage_histogram.png?raw=true" height="400">
<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_coverage_histogram_2.png?raw=true" height="400">
<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_coverage_histogram_3.png?raw=true" height="400">

<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_coverage_cumulative.png?raw=true" height="400">
<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_insert_size.png?raw=true" height="400">
<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/qualimap_gc_content.png?raw=true" height="400">

<img src="../output_WGBS/methylseq_V3_bwa_test/multiqc/bwameth/report_screenshots/picard_deduplication.png?raw=true" height="400">

## Analyzing methylation calls in methylkit

See 'code/09-MethylKit.Rmd'

## Analyzing CpGs

### note

If you don't have bedgraph outputs because you used --methylKit for running methylseq, this is the code to make them:

```
for bamFile in *.markdup.sorted.bam; do
        prefix=$(basename $bamFile .markdup.sorted.bam)

        #output bismark-style cytosine report
        singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --CHG --CHH --cytosine_report $reference $bamFile

        #output default bedgraph format (three, one for each context)
        singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --CHG --CHH $reference $bamFile

done
```

### Bedtools intersection of gene bodies (from gtf)

First, convert the gtf into a bed file:

```
salloc

cd references
module load gffread/0.12.7
module load bedtools2/2.31.1

gffread Pocillopora_acuta_HIv2.gtf --bed -o Pocillopora_acuta_HIv2.gtf.bed
awk '{OFS="\t"; print $1, $2, $3, $4}' Pocillopora_acuta_HIv2.gtf.bed > Pocillopora_acuta_HIv2.gtf.cleaned.bed.int
bedtools sort -i Pocillopora_acuta_HIv2.gtf.cleaned.bed.int > Pocillopora_acuta_HIv2.gtf.cleaned.bed

cd /home/zdellaert_uri_edu/zdellaert_uri_edu-shared/methylseq_bwa_fb/methyldackel
cp /home/zdellaert_uri_edu/LaserCoral/references/Pocillopora_acuta_HIv2.gtf.cleaned.bed .
```

Then: find the intersection between the CpG, CHG, and CHH sites and the gene bodies

```
nano scripts/methylation_efficiency.sh
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=scripts/outs_errs/"%x_error.%j"
#SBATCH --output=scripts/outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/

# load modules
module load bedtools2/2.31.1

# reference file
Pacuta_gff="Pocillopora_acuta_HIv2.gtf.cleaned.bed"

cd output_WGBS/methylseq_V3_bwa_test/methyldackel/

# loop through all samples - there is only one cytosine_report.txt file per sample

for file in *cytosine_report.txt; do

    sample=$(basename "$file" .markdup.sorted.cytosine_report.txt)

    # set output directory, sample-specific
    output_dir="${sample}_quant"
    mkdir -p $output_dir

    # name variables
    CpG_bed="${sample}.markdup.sorted_CpG.bedGraph"
    CHG_bed="${sample}.markdup.sorted_CHG.bedGraph"
    CHH_bed="${sample}.markdup.sorted_CHH.bedGraph"
    cytosine_report="${sample}.markdup.sorted.cytosine_report.txt"

    # extract CpG, CHH, and CHG methylation for each gene
    bedtools map -a $Pacuta_gff -b $CpG_bed -c 4 -o mean > $output_dir/gene_body_CpG_methylation.txt
    bedtools map -a $Pacuta_gff -b $CHG_bed -c 4 -o mean > $output_dir/gene_body_CHG_methylation.txt
    bedtools map -a $Pacuta_gff -b $CHH_bed -c 4 -o mean > $output_dir/gene_body_CHH_methylation.txt

    # extract total cytosines from cytosine_report
    awk '$6 == "CHG"' $cytosine_report > $output_dir/CHG_total.txt
    awk '$6 == "CHH"' $cytosine_report > $output_dir/CHH_total.txt

    # extract unmethylated cytosines from cytosine_report
    awk '$6 == "CHG" && $5 == 0' $cytosine_report > $output_dir/CHG_unmethylated.txt
    awk '$6 == "CHH" && $5 == 0' $cytosine_report > $output_dir/CHH_unmethylated.txt

done

### Coverage > 5x only

for file in *cytosine_report.txt; do

    sample=$(basename "$file" .markdup.sorted.cytosine_report.txt)

    # set output directory, sample-specific
    output_dir="${sample}_5Xquant"
    mkdir -p $output_dir

    # name variables
    CpG_bed="${sample}.markdup.sorted_CpG.bedGraph"
    CHG_bed="${sample}.markdup.sorted_CHG.bedGraph"
    CHH_bed="${sample}.markdup.sorted_CHH.bedGraph"
    cytosine_report="${sample}.markdup.sorted.cytosine_report.txt"

    # first calculate coverage by adding together methylated and unmethylated columns

    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5 + $6}' $CpG_bed > $output_dir/CpG_coverage.bed
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5 + $6}' $CHG_bed > $output_dir/CHG_coverage.bed
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5 + $6}' $CHH_bed > $output_dir/CHH_coverage.bed

    # filter to keep sites with coverage greater than 5x
    awk '$5 >= 5' $output_dir/CpG_coverage.bed > $output_dir/CpG_coverage_5x.bed
    awk '$5 >= 5' $output_dir/CHG_coverage.bed > $output_dir/CHG_coverage_5x.bed
    awk '$5 >= 5' $output_dir/CHH_coverage.bed > $output_dir/CHH_coverage_5x.bed

    # extract CpG, CHH, and CHG methylation for each gene
    bedtools map -a $Pacuta_gff -b $output_dir/CpG_coverage_5x.bed -c 4 -o mean > $output_dir/gene_body_CpG_methylation.txt
    bedtools map -a $Pacuta_gff -b $output_dir/CHG_coverage_5x.bed -c 4 -o mean > $output_dir/gene_body_CHG_methylation.txt
    bedtools map -a $Pacuta_gff -b $output_dir/CHH_coverage_5x.bed -c 4 -o mean > $output_dir/gene_body_CHH_methylation.txt

    # extract total cytosines from cytosine_report with coverage > 5x (calculate coverage by adding together methylated and unmethylated columns)
    awk '$6 == "CHG" && ($4 + $5) > 5' $cytosine_report > $output_dir/CHG_total_5x.txt
    awk '$6 == "CHH" && ($4 + $5) > 5' $cytosine_report > $output_dir/CHH_total_5x.txt

    # extract unmethylated cytosines from cytosine_report with coverage > 5x (calculate coverage by adding together methylated and unmethylated columns)
    awk '$6 == "CHG" && $5 == 0 && ($4 + $5) > 5' $cytosine_report > $output_dir/CHG_unmethylated_5x.txt
    awk '$6 == "CHH" && $5 == 0 && ($4 + $5) > 5' $cytosine_report > $output_dir/CHH_unmethylated_5x.txt

done
```

## Testing what happens if I filter for min-conversion efficeincy during the methylation extraction.

https://github.com/dpryan79/MethylDackel:

"Excluding partially converted reads
Some users wish to filter out reads that have evidence of incomplete bisulfite conversion. This can be determined by looking for CHH and CHG-context Cytosines in each read and determining their methylation state. Such filtering should generally be avoided, since there are regions in most genomes with at least some CHH and CHG-context Cytosine methylation, which would cause excess filtering in those regions. However, if you absolutely need to filter out only partially-converted alignments, you can do so with the --minConversionEfficiency option. The default is 0, meaning that CHH and CHG-context Cytosine conversion status is ignored. The maximum value is 1.0, meaning that 100% of the CHH and CHG-context Cytosines in an alignment must be converted to T. Note that an alignment with no CHH or CHG-context Cytosines will be given a conversion efficiency of 1.0."


This is recommended against, but I just want to see how it affects the calculations.

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=6
#SBATCH --mem=200GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"

module load samtools/1.19.2
module load apptainer/latest
reference=" ../../../references/Pocillopora_acuta_HIv2.assembly.fasta"

#singularity pull https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7

for bamFile in *.markdup.sorted.bam; do
	prefix=$(basename $bamFile .markdup.sorted.bam)

	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --CHG --CHH --cytosine_report --minConversionEfficiency 0.9 -o min_efficiency_test_new/min_90_${prefix} $reference $bamFile
	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --CHG --CHH --minConversionEfficiency 0.9 -o min_efficiency_test_new/min_90_${prefix} $reference $bamFile
	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --methylKit --minConversionEfficiency 0.9 -o min_efficiency_test_new/min_90_${prefix} $reference $bamFile
done
```

```
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH -p cpu -t 2:00:00
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"

# load modules
module load bedtools2/2.31.1

# reference file
Pacuta_gff="../../../../references/Pocillopora_acuta_HIv2.gtf.cleaned.bed"

# loop through all samples - there is only one cytosine_report.txt file per sample

cd min_efficiency_test_new

for file in *cytosine_report.txt; do

    sample=$(basename "$file" .cytosine_report.txt)

    # set output directory, sample-specific
    output_dir="${sample}_quant"
    mkdir -p $output_dir

    # name variables
    CpG_bed="${sample}_CpG.bedGraph"
    CHG_bed="${sample}_CHG.bedGraph"
    CHH_bed="${sample}_CHH.bedGraph"
    cytosine_report="${sample}.cytosine_report.txt"

    # extract CpG, CHH, and CHG methylation for each gene
    bedtools map -a $Pacuta_gff -b $CpG_bed -c 4 -o mean > $output_dir/gene_body_CpG_methylation.txt
    bedtools map -a $Pacuta_gff -b $CHG_bed -c 4 -o mean > $output_dir/gene_body_CHG_methylation.txt
    bedtools map -a $Pacuta_gff -b $CHH_bed -c 4 -o mean > $output_dir/gene_body_CHH_methylation.txt

    # extract conversion efficiency info

    #count unmethylated CHG
    awk '$4 == 0' $CHG_bed | wc -l > $output_dir/efficiency.txt

    #count unmethylated CHH
    awk '$4 == 0' $CHH_bed | wc -l >> $output_dir/efficiency.txt

    #count total CHG
    cat $CHG_bed | wc -l >> $output_dir/efficiency.txt

    #count total CHG
    cat $CHH_bed | wc -l >> $output_dir/efficiency.txt

    # will be calculated in R, but conversion efficiency = (unmethylated CHG + unmethylated CHH)/(total CHG + total CHH)
    # conversion efficiency = (efficiency.txt[1] + efficiency.txt[2])/(efficiency.txt[3] + efficiency.txt[4])
done
```

## Summary / Conclusions of where we stand with this data.

First: See [This markdown](https://github.com/zdellaert/LaserCoral/blob/main/code/10-CpG-Methylation-V3.md) for conversion efficiency info for the whole dataset, the dataset with 5X coverage, and the dataset filtered for > 90% conversion efficiency.

![alt text](10-CpG-Methylation-V3_files/figure-gfm/unnamed-chunk-11-1.png) ![alt text](10-CpG-Methylation-V3_files/figure-gfm/unnamed-chunk-11-2.png) ![alt text](10-CpG-Methylation-V3_files/figure-gfm/unnamed-chunk-11-3.png)

Two main files of DML/DMG interpretation:

1. [Whole Dataset, Mapped with BWA Meth (V3 trimming)](https://github.com/zdellaert/LaserCoral/blob/main/code/09-MethylKit.md)

Takeaways:

- There are 59 DMLs, 42 of which occur in gene bodies (DMGs)
- 4 DMGs are also DEGs
    - Pocillopora_acuta_HIv2___TS.g29592.t2
    - Pocillopora_acuta_HIv2___TS.g3401.t1
    - Pocillopora_acuta_HIv2___RNAseq.g1208.t1
    - Pocillopora_acuta_HIv2___TS.g19623.t2

| Query                                     | BaseMean  | Log2FoldChange | lfcSE  | PValue       | Padj        | QValue     | MethDiff  |
|-------------------------------------------|-----------|---------------|--------|-------------|------------|------------|-----------|
| Pocillopora_acuta_HIv2___TS.g3401.t1      | 92.69481  | -11.2018977   | 3.6354 | 6.46e-08    | 1.47e-06   | 0.02159195 | 15.20642  |
| Pocillopora_acuta_HIv2___TS.g19623.t2     | 269.40223 | -0.7395856    | 1.4279 | 1.10e-06    | 1.71e-05   | 0.04578924 | 15.57377  |
| Pocillopora_acuta_HIv2___RNAseq.g1208.t1  | 544.93022 | 0.7332125     | 1.4577 | 4.65e-05    | 4.29e-04   | 0.04856409 | 19.30342  |
| Pocillopora_acuta_HIv2___TS.g29592.t2     | 198.58417 | -0.4754287    | 1.1099 | 1.41e-03    | 7.88e-03   | 0.01725888 | 27.85388  |

2. [Same as above, but reads filtered for > 90% conversion efficiency](https://github.com/zdellaert/LaserCoral/blob/main/code/09-MethylKit-ConvFilt.md)

Takeaways:
- No DMGs exist in this dataset post-filtering, so none of these overlap with DEGs