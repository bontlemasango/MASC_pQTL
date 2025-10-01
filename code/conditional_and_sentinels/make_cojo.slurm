#!/bin/bash
#SBATCH --job-name=cojo
#SBATCH --output=logs/cojo_%A_%a.out
#SBATCH --error=logs/cojo_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

module load plink
module load gcta

MA_LIST=significant_sentinels.txt   # list of .ma files (one per line)
BFILE=/spaces/bmasango/bmasango/MASC_genotype/genoz
GCTA=/opt/exp_soft/bioinf/gcta/gcta-1.94.4

MA_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $MA_LIST)
base=$(basename "$MA_FILE" .ma)
IFS='_' read -r _ chr pos rest <<< "$base"

snplist="snps_${base}.txt"
out_plink="plink_subset_${base}"
out_cojo="cojo_${base}"

awk 'NR>1{print $1}' "$MA_FILE" > "$snplist"

plink --bfile $BFILE --extract "$snplist" --make-bed --out "$out_plink"

$GCTA --bfile "$out_plink" --chr "$chr" --cojo-file "$MA_FILE" --cojo-slct --out "$out_cojo"
rm -f ${out_plink}.bed ${out_plink}.bim ${out_plink}.fam
