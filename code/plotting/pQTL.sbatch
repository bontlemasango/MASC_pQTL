#!/bin/bash

#SBATCH --job-name=pQTL_plot
##SBATCH --account=gen1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
##SBATCH --array=1-37
#SBATCH --time=4:00:00
#SBATCH --output=slurm_out/%x-%A_%a.out
#SBATCH --export=None

eval "$(/home/bmasango/.local/bin/micromamba  shell hook --shell bash)"; micromamba activate regenie
python pQTL.py ../../data/processed/regenie_step2 hgnc_ucsc_b38_genes.tsv --genename
