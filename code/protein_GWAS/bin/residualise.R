#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
})

parser <- arg_parser(paste(
    'Residualise a Regenie-style phenotype(s) tab-delimited file',
    'on a separate covariates tab-delimited file.',
    'Files should be tab-delimited, with two starting columns: FID and IID.',
    'Headers are necessary.',
    'Each phenotype column will be residualised on ALL the given covariates.',
    '', sep='\n'))
parser %<>% add_argument('pheno', nargs=1, help='phenotypes tsv')
parser %<>% add_argument('covars', nargs=1, help='covariates tsv')
parser %<>% add_argument('--output', nargs=1, help='output file name', default='pheno.resid.tsv')
parser %<>% parse_args

cat(paste(
    'Running script with the following arguments',
    paste('  - pheno:', parser$pheno),
    paste('  - covars:', parser$covars),
    paste('  - output:', parser$output),
    '', sep='\n')
)

# # TESTING
# parser <- list(pheno='phenos.tsv', covars='covars.tsv', output='pheno.resid.tsv')
# # TESTING

# +
# Read the data
pheno <- read.table(parser$pheno, sep='\t', header=T)
covar <- read.table(parser$covars, sep='\t', header=T)

id_cols <- c('FID', 'IID')

# Convert FID/IID columns to uppercase in both files to ensure matching
names(pheno)[1:2] %<>% toupper
names(covar)[1:2] %<>% toupper

# Get the non-FID or IID columns
pheno_cols <- names(pheno)[3:ncol(pheno)]
covar_cols <- names(covar)[3:ncol(covar)]

# +
# Merge together to ensure rows are aligned
merged <- pheno %>% merge(covar,
                          all.x=T, sort=F,
                          by.x=id_cols,
                          by.y=id_cols)

# Get the residuals for each
for(phenotype in pheno_cols){
    lm_formula <- paste(covar_cols, collapse='+') %>%
                    paste0(phenotype, '~', .) %>%
                    as.formula
    phenotype_residuals <- lm_formula %>%
                                lm(data=merged, na.action='na.exclude') %>%
                                residuals
    merged[phenotype] <- phenotype_residuals
}

# Write out
out <- merged[, !(names(merged) %in% covar_cols)]
# -

write.table(out, file=parser$output, sep='\t', row.names=F, quote=F)
