# Trim Saga

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

1. Sliding window trimming
   1. *qtrim:                 WIN*
   2. **qtrim-win-size:        5**
2. Add right and lefthand polyG trimming
   1. -htrim-right: G -htrim-left: G

----

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
