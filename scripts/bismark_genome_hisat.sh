#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=100GB
#SBATCH -t 2:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load HISAT2/2.2.1-gompi-2022a

### Prepare Genome

mkdir -p references/Bisulfite_Genome_hisat2
cd references/Bisulfite_Genome_hisat2
cp ../Pocillopora_acuta_HIv2.assembly.fasta .

bismark_genome_preparation --hisat2 --verbose --parallel 10 ./
