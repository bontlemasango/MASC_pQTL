#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(readr)

#------------------------------
# INPUT FILES AND PARAMETERS
#------------------------------

sentinels_file <- "/spaces/bmasango/bmasango/proteomic_GWAS/data/processed/sentinels/sentinels_16_08_2025/significant_sentinels.txt"
regenie_dir <- "."

# Output directory
out_dir <- "../all_ma_files/"
if (!dir.exists(out_dir)) dir.create(out_dir)

window <- 1e6  # +/-1Mb

#------------------------------
# GET ARRAY INDEX FROM SLURM
#------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: Must provide sentinel index as argument")
}
i <- as.numeric(args[1])

#------------------------------
# LOAD SENTINELS
#------------------------------
sentinels <- read_tsv(sentinels_file, show_col_types = FALSE)

if (i > nrow(sentinels)) stop("Index exceeds number of sentinels")

#------------------------------
# GET INFO FOR THIS SENTINEL
#------------------------------
chr <- sentinels$chrom[i]
pos <- sentinels$pos[i]
rsid <- sentinels$rsid[i]
protein <- sentinels$protein[i]

out_file <- file.path(out_dir, paste0("sentinel_", chr, "_", pos, "_", protein, ".ma"))

# Skip if already done
if (file.exists(out_file)) {
  cat("Skipping already completed sentinel:", rsid, "->", out_file, "\n")
  quit(status = 0)
}

#------------------------------
# FIND REGENIE FILE
#------------------------------
regenie_file <- list.files(
  path = regenie_dir,
  pattern = paste0("chr_", chr, ".*", protein, ".*\\.regenie$"),
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(regenie_file) == 0) {
  cat("WARNING: File not found for sentinel:", rsid, "chr", chr, protein, "\n")
  quit(status = 0)
}

regenie_file <- regenie_file[1]

#------------------------------
# READ ONLY NECESSARY COLUMNS
#------------------------------
gwas <- fread(
  regenie_file,
  select = c("ID","ALLELE0","ALLELE1","A1FREQ","BETA","SE","LOG10P","N")
)
gwas <- as.data.frame(gwas)
names(gwas) <- toupper(names(gwas))

# Convert GENPOS if exists
if ("GENPOS" %in% names(gwas)) {
  if (!is.numeric(gwas$GENPOS)) {
    gwas$GENPOS <- suppressWarnings(as.numeric(as.character(gwas$GENPOS)))
  }
} else {
  gwas$GENPOS <- NA
}

#------------------------------
# SUBSET +/-1Mb
#------------------------------
if (!is.na(gwas$GENPOS[1])) {
  idx <- which(gwas$GENPOS >= pos - window & gwas$GENPOS <= pos + window)
} else {
  idx <- seq_len(nrow(gwas))  # take all if no GENPOS
}

if (length(idx) == 0) {
  cat("No SNPs found for sentinel:", rsid, "at chr", chr, pos, "\n")
  quit(status = 0)
}

region <- gwas[idx, , drop = FALSE]

#------------------------------
# CONVERT LOG10P -> P AND FORMAT
#------------------------------
region <- region %>%
  mutate(P = 10^(-LOG10P)) %>%
  transmute(
    SNP = ID,
    A1 = ALLELE1,
    A2 = ALLELE0,
    freq = A1FREQ,
    b = BETA,
    se = SE,
    p = P,
    N = N
  )

#------------------------------
# WRITE OUTPUT
#------------------------------
write_tsv(region, out_file)
cat("Wrote:", out_file, "with", nrow(region), "SNPs\n")
