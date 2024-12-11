#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=500GB
#SBATCH --error="%x_error.%j" 
#SBATCH --output="%x_output.%j" 
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

SINGULARITY_IMAGE="docker://memesuite/memesuite:latest"

cd ../output_RNA/differential_expression/TFs

module load apptainer/latest

singularity exec --cleanenv --debug $SINGULARITY_IMAGE meme promoters_500_upstream.fasta -maxw 25 -mod anr -nmotifs 2 -oc meme_output
