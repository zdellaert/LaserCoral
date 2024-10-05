#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/output_RNA #set working directory

# load modules needed
module load StringTie/2.1.4-GCC-9.3.0 #stringtie module, includes Python 3.8.5

# move into stringtie directory
cd stringtie

# make input file
for filename in *.gtf; do
    sample_name=$(echo "$filename" | awk -F'[.]' '{print $1}')

    echo $sample_name $PWD/$filename
done > listGTF.txt

#Compile the gene count matrix
python ../../scripts/prepDE.py3 -g LCM_RNA_gene_count_matrix.csv -i listGTF.txt

echo "Gene count matrix compiled." $(date)
