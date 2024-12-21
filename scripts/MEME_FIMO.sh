#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
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

# run FIMO 

#singularity exec --cleanenv $SINGULARITY_IMAGE fimo --oc fimo_output_upOralEpi meme_output_upOralEpi/motifs.meme --motif "CCGCCATBTTK" promoters_500_upstream.fasta
singularity exec --cleanenv $SINGULARITY_IMAGE fimo --oc fimo_output_upOralEpi --motif "CCGCCATBTTK" meme_output_upOralEpi/meme.txt promoters_500_upstream.fasta
