#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=12
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=100GB
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

SINGULARITY_IMAGE="docker://memesuite/memesuite:latest"

cd ../output_RNA/differential_expression/TFs

module load apptainer/latest

singularity exec --cleanenv $SINGULARITY_IMAGE meme promoters_500_upstream_upAboral.fasta -dna -maxw 25 -mod anr -nmotifs 10 -p 6 -oc meme_output_upAboral -revcomp

singularity exec --cleanenv $SINGULARITY_IMAGE meme promoters_500_upstream_upOralEpi.fasta -dna -maxw 25 -mod anr -nmotifs 10 -p 6 -oc meme_output_upOralEpi -revcomp

