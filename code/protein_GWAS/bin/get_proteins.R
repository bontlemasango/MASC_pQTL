#!/usr/bin/env Rscript
# # Code to retrieve the EXCEED protein data and format it appropriately for REGENIE

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
    'Script to create genetic-data-compatible protein phenotype file',
    'Creates tab-delimited file, with two starting columns: FID and IID.',
    'Additionally creates a mapping file to convert the individual IDs from DIA-NN format to genetic format',
    '', sep='\n'))
parser %<>% add_argument('--output', nargs=1, help='output file name', default='proteins.tsv')
parser %<>% parse_args

proteins <- readxl::read_excel("/spaces/bmasango/bmasango/proteomic_GWAS/data/raw/masc_baseline_data.xlsx")
all_ids <- readxl::read_excel("/spaces/bmasango/bmasango/proteomic_GWAS/data/raw/Merged GSK,AWI-Gen,SWEET & H3A IDs_lisa_25June2020.xlsx")
all_ids <- all_ids %>% select(`AWI-Gen ID`,gsk_id)
#do this if you want to select more than 1 protein the : thing i.e (proteins <- proteins %>% select(gsk_id,`BMP-6`:MEPE))
proteins <- proteins %>% select(gsk_id,`BMP-6`:`MEPE`)
protein_GWAS <- merge(all_ids, proteins, by = "gsk_id") %>%
  mutate(FID=`AWI-Gen ID`, IID=`AWI-Gen ID`)  %>% 
 select(FID, IID, everything()) %>%
select(-`AWI-Gen ID`,-gsk_id)
write_tsv(protein_GWAS,parser$output)
