#! /bin/env Rscript
# Note: run using Alice R module (rtracklayer version )
# Collect all the UCSC data on gene/protein position
library(tidyverse)
library(rtracklayer)

get_ucsc_table <- function(file.out, track, tbl){
    if(file.exists(file.out)){
        print(paste(file.out, 'already exists'))
        return()
    }
    ucsc = browserSession("UCSC")
    genome(ucsc) <- "hg38"
    
    refseq <- ucsc %>%
                ucscTableQuery(track=track, table=tbl) %>%
                (rtracklayer::getTable) %>%
                as_tibble

    write.table(refseq, file=file.out, sep='\t')
}

# Pull data for both the 'hgnc' gene data and the ncbi gene data
get_ucsc_table('tmp_data/hgnc_ucsc_b38_genes.tsv', 'hgnc', 'hgnc')
get_ucsc_table('tmp_data/ncbi_ucsc_b38_genes.tsv', 'NCBI RefSeq', 'refGene')