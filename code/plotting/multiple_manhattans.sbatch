#!/bin/bash

#SBATCH --job-name=manhattans
##SBATCH --account=gen1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-37
#SBATCH --time=1:00:00
#SBATCH --output=slurm_out/%x-%A_%a.out
#SBATCH --export=None


eval "$(/home/bmasango/.local/bin/micromamba  shell hook --shell bash)"; micromamba activate regenie
python multiple_manhattans.py ../../data/processed/regenie_step2 ${SLURM_ARRAY_TASK_ID} 
