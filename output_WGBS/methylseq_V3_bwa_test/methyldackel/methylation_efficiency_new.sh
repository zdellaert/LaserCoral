#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --qos=short
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
