#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=2
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH --error=../scripts/outs_errs/"%x_error.%j"
#SBATCH --output=../scripts/outs_errs/"%x_output.%j"
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

SINGULARITY_IMAGE="docker://memesuite/memesuite:latest"

#Download Motif databases
#cd ../references
#wget https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.25.tgz
#tar -xvzf motif_databases.12.25.tgz
#mv motif_databases motif_dbs

cd ../output_RNA/differential_expression/TFs

module load apptainer/latest

# Run STREME for relative enrichment, only report motifs with e-value < 0.05
singularity exec --cleanenv $SINGULARITY_IMAGE streme --oc streme_output_upAboral \
  --p promoters_500_upstream_upAboral.fasta \
  --n promoters_500_upstream.fasta \
  --patience 50 \
  --dna --thresh 0.05

singularity exec --cleanenv $SINGULARITY_IMAGE streme --oc streme_output_upOralEpi \
  --p promoters_500_upstream_upOralEpi.fasta \
  --n promoters_500_upstream.fasta \
  --patience 50 \
  --dna --thresh 0.05 

# run TOMTOM on the MEME-identified motifs from above

singularity exec --cleanenv $SINGULARITY_IMAGE tomtom -no-ssc -oc tomtom_streme_output_upAboral -min-overlap 5 -dist pearson -thresh 0.05 streme_output_upAboral/streme.txt ../../../references/motif_dbs/JASPAR/JASPAR2022_CORE_non-redundant_v2.meme ../../../references/motif_dbs/EUKARYOTE/homeodomain.meme ../../../references/motif_dbs/EUKARYOTE/jolma2013.meme

singularity exec --cleanenv $SINGULARITY_IMAGE tomtom -no-ssc -oc tomtom_streme_output_upOralEpi -min-overlap 5 -dist pearson -thresh 0.05 streme_output_upOralEpi/streme.txt ../../../references/motif_dbs/JASPAR/JASPAR2022_CORE_non-redundant_v2.meme ../../../references/motif_dbs/EUKARYOTE/homeodomain.meme ../../../references/motif_dbs/EUKARYOTE/jolma2013.meme
