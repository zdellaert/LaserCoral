#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/references

# load modules needed
module load gffread/0.12.7-GCCcore-11.2.0

# "Clean GFF" file if necessary, then convert cleaned file into a GTF

gffread -E Pocillopora_acuta_HIv2.genes.gff3 -T -o Pocillopora_acuta_HIv2.gtf 

echo "Check how many fields are in each row of the file; currently there are rows with two different lenghts: 10 and 12"
awk '{print NF}' Pocillopora_acuta_HIv2.gtf | sort -u

# Use awk to add "gene_id = TRANSCRIPT_ID" to each of the rows that only have a transcript id listed (the non-transcript lines)
awk 'BEGIN {FS=OFS="\t"} {if ($9 ~ /transcript_id/ && $9 !~ /gene_id/) {match($9, /transcript_id "([^"]+)";/, a); $9 = $9 " gene_id \"" a[1] "\";"}; print}' Pocillopora_acuta_HIv2.gtf > Pocillopora_acuta_HIv2_modified.gtf

echo "Check how many fields are in each row of the file; Now all the rows should be the same length and only one value should be printed, 12"
awk '{print NF}' Pocillopora_acuta_HIv2_modified.gtf | sort -u

# remove the non-modified file
rm Pocillopora_acuta_HIv2.gtf

# rename the modified file
mv Pocillopora_acuta_HIv2_modified.gtf Pocillopora_acuta_HIv2.gtf
