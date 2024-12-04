#!/bin/bash

#SBATCH --job-name=GWAS_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --export=NONE
#SBATCH --account=gen1

shopt -s expand_aliases
source ~/.bashrc

if [ "$1" == 'dag' ]; then
    nextflow \
    -log dry_run.log \
    run workflow.nf \
    -c workflow.config \
    -preview \
    -with-dag ../../results/figures/gwas.dag.svg

else
    nextflow \
    -log workflow.log \
    run workflow.nf \
    -c workflow.config \
    -with-trace \
    -resume
fi
