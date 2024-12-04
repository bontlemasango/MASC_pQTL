#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
})

parser <- arg_parser(paste(
    'Rank inverse-normal transform a Regenie-style phenotype(s) tab-delimited file',
    'File should be tab-delimited, with two starting columns: FID and IID.',
    'Headers are necessary.',
    'Each phenotype column will be transformed separately, NAs ignored.',
    '', sep='\n'))
parser %<>% add_argument('pheno', nargs=1, help='phenotypes tsv')
parser %<>% add_argument('--output', nargs=1, help='output file name', default='pheno.INT.tsv')
parser %<>% parse_args

cat(paste(
    'Running script with the following arguments',
    paste('  - pheno:', parser$pheno),
    paste('  - output:', parser$output),
    '', sep='\n')
)

# # TESTING
# parser <- list(pheno='pheno.resid.tsv', output='pheno.INT.tsv')
# # TESTING

# +
# Read the data
pheno <- read.table(parser$pheno, sep='\t', header=T)

id_cols <- c('FID', 'IID')

# Convert FID/IID columns to uppercase in file to ensure matching
names(pheno)[1:2] %<>% toupper

# Get the non-FID or IID columns
pheno_cols <- names(pheno)[3:ncol(pheno)]

# +
# Convert the values to ranks
out <- pheno
for(phenotype in pheno_cols){
    out[phenotype] <- rank(pheno[phenotype], ties.method='average', na.last=T)
    # Re-insert NAs where necessary
    out[phenotype][is.na(pheno[phenotype])] <- NA
}

# Convert the ranks to z-scores
for(phenotype in pheno_cols){
    x <- out[[phenotype]]
    # Inverse normal transform
    out[phenotype] <- qnorm((x-0.5)/max(x, na.rm=T))
}
# -

write.table(out, file=parser$output, sep='\t', row.names=F, quote=F)
