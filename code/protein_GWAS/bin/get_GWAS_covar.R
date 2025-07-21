#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
    library(dplyr)
    library(tidyverse)
})

parser <- arg_parser(paste(
    'Create the file to get the covariates (PCs and batches) for EXCEED GWAS.',
    'Takes 1 input:',
    '  output filename, default: exceed.eur.geno.covar',
    'Output is a tab-delimited files, with the first columns being FID and IID.',
    '', sep='\n'))
parser %<>% add_argument('--output', nargs=1, help='output file name', default='exceed.eur.geno.covar')
parser %<>% parse_args

AWI10pcs <- read_table("/spaces/bmasango/bmasango/proteomic_GWAS/data/raw/AWI10pcs.txt")
write_tsv(AWI10pcs,parser$output)
