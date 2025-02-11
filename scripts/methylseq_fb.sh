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

# nextflow run nf-core/methylseq \
#   --input ./data_WGBS/LCM_methylseq_input_fb.csv \
#   --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_bwa_test \
#   --aligner bwameth \
#   --methyl_kit \
#   --igenomes_ignore \
#   --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
#   --relax_mismatches \
#   --save_align_intermeds \
#   --skip_trimming \
#   --email zdellaert@uri.edu \
#   -profile unity \
#   -name methylseq_bwa_test

nextflow run nf-core/methylseq \
  --input ./data_WGBS/LCM_methylseq_input_fb.csv \
  --outdir /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_bwa_test \
  --aligner bwameth \
  --methyl_kit \
  --igenomes_ignore \
  --fasta ./references/Pocillopora_acuta_HIv2.assembly.fasta \
  --fasta_index ./references/Pocillopora_acuta_HIv2.assembly.fasta.fai \
  --bwameth_index ./references/BwamethIndex/ \
  --run_qualimap \
  --relax_mismatches \
  --save_align_intermeds \
  --skip_trimming \
  --email zdellaert@uri.edu \
  -profile unity \
  -resume 98b8ae7d-d97b-4342-ac97-f48ec63a4bcc
