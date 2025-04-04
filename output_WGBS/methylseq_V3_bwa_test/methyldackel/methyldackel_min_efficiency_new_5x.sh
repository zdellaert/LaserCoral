#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=6
#SBATCH --mem=200GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"

## This is an alternative script to run the bwa-meth post processing manually if methylseq quits after the alignment and doesn't want to resume....

module load samtools/1.19.2
module load apptainer/latest
reference=" ../../../references/Pocillopora_acuta_HIv2.assembly.fasta"

#singularity pull https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7

for bamFile in /scratch3/workspace/zdellaert_uri_edu-shared/methylseq_V3_bwa_test/bwameth/deduplicated/LCM_*.markdup.sorted.bam; do
	prefix=$(basename $bamFile .markdup.sorted.bam)

	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --CHG --CHH --minDepth 5 --cytosine_report --minConversionEfficiency 0.9 -o min_efficiency_test_new_5x/min_90_${prefix} $reference $bamFile
	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --CHG --CHH --minDepth 5 --minConversionEfficiency 0.9 -o min_efficiency_test_new_5x/min_90_${prefix} $reference $bamFile
	singularity exec methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 1 --methylKit --minDepth 5 --minConversionEfficiency 0.9 -o min_efficiency_test_new_5x/min_90_${prefix} $reference $bamFile
done
