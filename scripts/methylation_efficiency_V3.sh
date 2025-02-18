#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"
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
