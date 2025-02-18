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
