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
