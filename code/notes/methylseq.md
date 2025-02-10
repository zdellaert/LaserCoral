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


## Run with bwa-meth pipeline

Start with flexbar-trimmed data.

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
#SBATCH --nodes=1 --ntasks-per-node=30
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=500GB
#SBATCH -t 48:00:00
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



## on V3...

Start with flexbar-trimmed data.

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
   --relax_mismatches \
   --save_align_intermeds \
   --skip_trimming \
   --email zdellaert@uri.edu \
   -profile unity \
   -name methylseq_V3_bwa_test
```

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
    
    if [ ! -f "/scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/alignments/${prefix}.bam" ]; then
    echo "ERROR: Expected BAM file '/scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/alignments/${prefix}.bam' not found!" >&2
    exit 1
    fi

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


## Methylseq BWA-post processing script (if methylseq doesn't want to work)

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

	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract $reference --methylKit $bamFile.dedup.bam

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