#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=48
#SBATCH --mem=300GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

## This is an alternative script to run the bwa-meth post processing manually if methylseq quits after the alignment and doesn't want to resume....

module load samtools/1.19.2
module load apptainer/latest
reference="Pocillopora_acuta_HIv2.assembly.fasta"

cd output_WGBS/methylseq_V3_bwa_test/bwameth/deduplicated/

for bamFile in *.markdup.sorted.bam; do
	prefix=$(basename $bamFile .markdup.sorted.bam)

	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --CHG --CHH --cytosine_report $reference $bamFile
	singularity exec https://depot.galaxyproject.org/singularity/methyldackel:0.6.1--he4a0461_7 MethylDackel extract -@ 24 --CHG --CHH $reference $bamFile
done

mv *bedGraph ../../methyldackel
mv *cytosine_report.txt ../../methyldackel

