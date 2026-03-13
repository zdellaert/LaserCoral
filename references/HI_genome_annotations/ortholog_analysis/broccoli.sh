#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=20 #split one task over multiple CPU
#SBATCH --mem=500GB
#SBATCH -p cpu-preempt
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file

#load conda and activate conda environment
module load conda/latest
conda activate /work/pi_hputnam_uri_edu/conda/envs/env-broccoli

#load additional programs needed to run broccoli
module load uri/main
module load diamond/2.1.7
module load all/FastTree/2.1.11-GCCcore-12.3.0

mkdir -p output/broccoli

cd output/broccoli
python ~/broccoli.py -dir ../../cnidarian_protein_fastas/ -ext '.faa' -path_fasttree FastTree -threads 8
