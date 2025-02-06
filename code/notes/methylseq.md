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
