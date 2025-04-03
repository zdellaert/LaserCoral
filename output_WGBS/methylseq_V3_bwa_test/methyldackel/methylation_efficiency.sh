#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --mem=16G --qos=short
#SBATCH -p cpu -t 2:00:00
#SBATCH --error=outs_errs/"%x_error.%j"
#SBATCH --output=outs_errs/"%x_output.%j"

# load modules
module load bedtools2/2.31.1

# reference file
Pacuta_gff="../../../references/Pocillopora_acuta_HIv2.gtf.cleaned.bed"

# loop through all samples - there is only one cytosine_report.txt file per sample

## 5X coverage

for file in *cytosine_report.txt; do

    sample=$(basename "$file" .markdup.sorted.cytosine_report.txt)

    # set output directory, sample-specific
    output_dir="${sample}_5Xquant"

    # name variables
    CpG_bed="${output_dir}/CpG_coverage_5x.bed"
    CHG_bed="${output_dir}/CHG_coverage_5x.bed"
    CHH_bed="${output_dir}/CHH_coverage_5x.bed"

    # extract conversion efficiency info
    awk '$4 == 0' $CHG_bed | wc -l > $output_dir/efficiency.txt
    awk '$4 == 0' $CHH_bed | wc -l >> $output_dir/efficiency.txt
    cat $CHG_bed | wc -l >> $output_dir/efficiency.txt
    cat $CHH_bed | wc -l >> $output_dir/efficiency.txt
done
