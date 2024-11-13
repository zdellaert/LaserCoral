#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500GB
#SBATCH --export=NONE
#SBATCH --error=../scripts/outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=../scripts/outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH -D /project/pi_hputnam_uri_edu/zdellaert/LaserCoral/references #set working directory

## Warning for putnam lab users of this script, this is written to run on Unity, not andromeda

# Use the gtf file run in Stringtie above:

# the gtf has to have the gene and transcript IDs be different, so append -T to each transcript_id and save to a new file
sed 's/transcript_id "\([^"]*\)"/transcript_id "\1-T"/g' Pocillopora_acuta_HIv2.gtf > Pocillopora_acuta_HIv2_modified.gtf

# activate environment
module load miniconda/22.11.1-1 #load miniconda
conda activate /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext

#mkdir -p geneext
cd geneext

# use --clip_strand both to not allow GeneExt to create overlaps on the same strand

python /work/pi_hputnam_uri_edu/conda/envs/GeneExt/geneext.py --verbose=3 \
    -g ../Pocillopora_acuta_HIv2_modified.gtf \
    -b ../../output_RNA/hisat2/merge.bam \
    -o Pocillopora_acuta_GeneExt.gtf \
    -j $SLURM_CPUS_PER_TASK --clip_strand both --force
