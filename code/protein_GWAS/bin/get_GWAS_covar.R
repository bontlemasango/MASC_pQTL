#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
    library(dplyr)
})

parser <- arg_parser(paste(
    'Create the file to get the covariates (PCs and batches) for EXCEED GWAS.',
    'Takes 1 input:',
    '  output filename, default: exceed.eur.geno.covar',
    'Output is a tab-delimited files, with the first columns being FID and IID.',
    '', sep='\n'))
parser %<>% add_argument('--output', nargs=1, help='output file name', default='exceed.eur.geno.covar')
parser %<>% parse_args

cat(paste(
    'Running script with the following arguments',
    paste('  - output:', parser$output),
    '', sep='\n')
)

# +
## Get PCs
pcs.eur <- '/rfs/EXCEED/Genotyping/freeze2/ancestry/EUR_PCA.eigenvec'

pcs.eur %<>% read.table(sep=' ', header=F) %>%
                set_names(c('FID', 'IID', paste0('PC', 1:(ncol(.)-2))))
# -


## Get the batch (freeze 1 and freeze 2)
pcs.eur$batch <- gsub('.*_(?=b)', '', pcs.eur$IID, perl=T)

# Some Europeans have duplicate genotypes (b0dup), we will remove these
pcs.eur %<>% subset(batch!='b0dup')

# Write out the covariate files out
write.table(pcs.eur, file=parser$output, sep='\t', row.names=F, quote=F)
print(paste('Written European PCs and batches to:', parser$output))
