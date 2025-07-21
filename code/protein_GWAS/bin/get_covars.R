#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
    library(dplyr)
    library(stringr)
    library(readxl)
    library(tidyverse)
})

# +

parser <- arg_parser(paste(
    'Script to get the covariates for the protein phenotype file',
    'Creates tab-delimited file, with two starting columns: FID and IID.',
    'Additionally creates a mapping file to convert the individual IDs from DIA-NN format to genetic format',
    '', sep='\n'))
parser %<>% add_argument('--output', nargs=1, help='output file name', default='covars.tsv')
parser %<>% parse_args()

proteins <- readxl::read_excel("/spaces/bmasango/bmasango/proteomic_GWAS/data/raw/masc_baseline_data.xlsx")
all_ids <- readxl::read_excel("/spaces/bmasango/bmasango/proteomic_GWAS/data/raw/Merged GSK,AWI-Gen,SWEET & H3A IDs_lisa_25June2020.xlsx")
all_ids <- all_ids %>% select(`AWI-Gen ID`,gsk_id)
covariates_protein <- proteins %>% select(gsk_id,age,gender)
covariates_protein <- merge(all_ids, covariates_protein, by = "gsk_id")
covariates_protein <- covariates_protein %>% mutate(FID=`AWI-Gen ID`, IID=`AWI-Gen ID`)
covariates_protein <- covariates_protein %>% select(FID, IID, everything())
covariates_protein <- covariates_protein %>% select(-`AWI-Gen ID`,-gsk_id)
write_tsv(covariates_protein,parser$output)
