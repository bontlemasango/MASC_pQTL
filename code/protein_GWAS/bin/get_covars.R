#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
    library(dplyr)
    library(stringr)
    library(readxl)
})

# +
# Necessary files
sex_dob <- '/rfs/EXCEED/PrimaryCare/temporary_for_epistat/arranged_sex_dob_21-02-2024.csv'
ms_batch <- '/rfs/EXCEED/Proteomics/20230216_sample_metadata.xlsx'

parser <- arg_parser(paste(
    'Script to get the covariates for the protein phenotype file',
    'Creates tab-delimited file, with two starting columns: FID and IID.',
    'Additionally creates a mapping file to convert the individual IDs from DIA-NN format to genetic format',
    '', sep='\n'))
parser %<>% add_argument('pheno', nargs=1, help='input phenotype file, with first columsn FID, IID. Tab-delimited')
parser %<>% add_argument('mapping-file', nargs=1, help='Name of mapping file to read. Tab-delimited')
parser %<>% add_argument('--output', nargs=1, help='output file name', default='covars.tsv')
parser %<>% parse_args()

cat(paste(
    'Running script with the following arguments',
    paste('  - pheno:', parser$pheno),
    paste('  - mapping-file:', parser$mapping_file),
    paste('  - output:', parser$output),
    '', sep='\n')
)
# -

sex_dob %<>% read.csv %>% select(exceed_id, sex, dob)
ms_batch %<>% read_excel %>% select(File.Name, MS.Batch, MS.Date)
pheno <- read.table(parser$pheno, sep='\t', header=T) %>% select(FID, IID)
mapping_file <- read.table(parser$mapping_file, sep='\t', header=T) %>% unique

# +
covars <- pheno %>%
            merge(mapping_file, by.x='IID', by.y='plink_id', all.x=T) %>%
            merge(ms_batch, by.x='diann_id', by.y='File.Name', all.x=T) %>%
            merge(sex_dob, by='exceed_id', all.x=T)

covars$age <- as.Date(covars$MS.Date, format='%Y%m%d') %>%
                subtract(as.Date(covars$dob, format='%Y-%m-%d')) %>%
                as.numeric %>%
                divide_by_int(365.25)

# Check the covariates file is unique individuals
covars$IID %>% duplicated %>% sum %>% paste('Duplicated individuals:',.) %>% print 

# +
# Remove people younger than 18 - incorrect dob
covars$age[which(covars$age<18)] <- NA

# Map the MS batches to alphabetical (to ensure it is considered unordered discrete)
covars$batch <- sapply(covars$MS.Batch, FUN=function(x) intToUtf8(x + 96))

# Select out the columns we want
covars %<>% select(FID, IID, sex, age, batch)
# -

write.table(covars, file=parser$output, row.names=F, quote=F, sep='\t')
