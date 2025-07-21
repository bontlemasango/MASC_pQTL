#!/bin/bash
#SBATCH --job-name=extract_log10p
#SBATCH --output=extract_log10p.out
#SBATCH --error=extract_log10p.err
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Load any necessary modules 
# module load python/3.9

# Activate Python environment
# source ~/myenv/bin/activate

# Run the script with arguments
python /spaces/bmasango/bmasango/proteomic_GWAS/code/plotting/extract_minlog10p.py /spaces/bmasango/bmasango/proteomic_GWAS/data/processed/regenie_step2

