#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=50GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=zdellaert@uri.edu #your email to send notifications
#SBATCH -D /data/putnamlab/zdellaert/LaserCoral/references #set working directory

# Use the gtf file run in Stringtie above:

# the gtf has to have the gene and transcript IDs be different, so append -T to each transcript_id and save to a new file
sed 's/transcript_id "\([^"]*\)"/transcript_id "\1-T"/g' Pocillopora_acuta_HIv2.gtf > Pocillopora_acuta_HIv2_modified.gtf

# activate environment
source ~/.bashrc
module load Miniconda3/4.9.2
conda activate geneext

mkdir -p geneext
cd geneext

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand

python /data/putnamlab/conda/GeneExt/geneext.py --verbose=3 \
    -g ../Pocillopora_acuta_HIv2_modified.gtf \
    -b ../../output_RNA/hisat2/merge.bam \
    -o Pocillopora_acuta_GeneExt.gtf \
    -j 18 --clip_strand both
