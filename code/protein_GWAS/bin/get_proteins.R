#!/usr/bin/env Rscript
# # Code to retrieve the EXCEED protein data and format it appropriately for REGENIE

suppressPackageStartupMessages({
    library(magrittr)
    library(argparser)
    library(dplyr)
    library(stringr)
    library(readxl)
    library(exceedapi)
})

# Necessary files
proteins <- '/rfs/EXCEED/Proteomics/DiannResults-ProjectSpecificLibrary/report.pg_matrix.tsv'
to_labid <- '/rfs/EXCEED/Proteomics/20230216_sample_metadata.xlsx'
ancestry <- '/rfs/EXCEED/Genotyping/freeze2/ancestry/EXCEED_freeze2_clustered_ancestry_kmeans.txt'
plink.fam <- '/scratch/gen1/bl212/exceed/proteomic_GWAS_freeze2/data/raw/tmp/EXCEED_freeze2_genotyped_TOPMed_b38.fam'

# +
parser <- arg_parser(paste(
    'Script to create genetic-data-compatible protein phenotype file',
    'Creates tab-delimited file, with two starting columns: FID and IID.',
    'Additionally creates a mapping file to convert the individual IDs from DIA-NN format to genetic format',
    '', sep='\n'))
parser %<>% add_argument('--output', nargs=1, help='output file name', default='proteins.tsv')
parser %<>% add_argument('--mapping-file', nargs=1, help='Name of mapping file to read. If it does not exist, one will be created.', default='protein_map.tsv')
parser %<>% parse_args

cat(paste(
    'Running script with the following arguments',
    paste('  - output:', parser$output),
    paste('  - mapping-file:', parser$mapping_file),
    '', sep='\n')
)
# -

# ## Step 1 Read the proteins in and format so each row is an individual and each column is a protein

# +
# Read the proteins dataframe
proteins %<>% read.csv(sep='\t', check.names=FALSE)

# Store the new column names and individual IDS
proteins.colnames <- proteins %>% pull(Protein.Group) %>% gsub(';','_', .)
proteins.ids <- names(proteins) %>% subset(grepl('/s-', .))

# Transpose the dataframe so proteins are columns and each row is a different individual
proteins %<>% select(all_of(proteins.ids)) %>%
                t %>%
                data.frame(row.names=NULL) %>%
                set_names(proteins.colnames) %>%
                mutate(diann_id=proteins.ids)
# -

# ## Step 2 Create a mapping file to convert the DIA-NN IDs (individual IDs) to genomic barcode IDs 

if(file.exists(parser$mapping_file)){
    mapping_file <- read.csv(parser$mapping_file, sep='\t')
}else{
    # Read the DIA-NN ID to Lab ID dataframe
    to_labid %<>% read_excel() %>% select(diann_id=File.Name, lab_id=`Sample ID`)

    # EXCEEDapi stores LabID-to-ThrivaID, ThrivaID-to-PsuedoID,
    # PsuedoID-to-UUID, UUID-to-ExceedID, ExceedID-to-SalivaBarcode
    exceed <- exceedapi::exceed_client(); exceed$auth()

    # Collect all the required mapping datasets
    to_thrivaid <- exceed$datastore('thriva', 'antibodytest-returned-samples') %>%
                    collect %>%
                    select(lab_id=Laboratory_ID, thriva_id=THR_Number)
    to_pid <- NULL
    for(file_id in 2:4){
        tmp <- exceed$datastore('thriva', 'antibodytest-specimen-ids', file_id) %>% collect
        to_pid %<>% rbind(tmp)
    }
    to_pid %<>% select(thriva_id=Tests_Sample_ID, pid=Tests_External_User_ID)
    to_uuid <- exceed$identities('thriva') %>% 
                collect %>%
                select(pid=pid, uuid=uuid)
    to_exceedid <- rbind(collect(exceed$identities('exceed')),
                     collect(exceed$identities('exceed_civicrm'))) %>%
                    select(uuid=uuid, exceed_id=pid) %>%
                    mutate(exceed_id=str_pad(exceed_id, 6, pad=0)) # Add padding if necessary
    to_barcode <- exceed$specimens() %>%
                    collect %>%
                    select(exceed_id=exceed_id, specimen_id=specimen_id) %>%
                    mutate(exceed_id=str_pad(exceed_id, 6, pad=0))

    # Convert the saliva barcodes to the IDs used in the genomic data
    plink.fam %<>% read.csv(sep=' ', header=F) %>%
                    select(V2) %>%
                    set_names('plink_id') %>%
                    mutate(specimen_id=gsub('_b.*', '', plink_id))

    # Create the necessary mapping file
    mapping_file <- to_labid %>%
                    merge(to_thrivaid, all.x=T) %>%
                    merge(to_pid, all.x=T) %>%
                    merge(to_uuid, all.x=T) %>%
                    merge(to_exceedid, all.x=T) %>%
                    merge(to_barcode, all.x=T) %>%
                    merge(plink.fam, all.x=T)
    # Re-order the columns so they are in the order they've been mapped
    mapping_file %<>% select(diann_id, thriva_id, pid, uuid, exceed_id, specimen_id, plink_id)
}

# +
# Note: our mapping file has many duplicated rows.
# This is due to uuids mapping to multiple Exceed IDs
print(paste('Number of rows of mapping file:', nrow(mapping_file)))

# We will check that the individual IDs from the protein dataset are
# a 1-to-1 mapping for the genetic data
diannid_to_plinkid <- mapping_file %>% 
                        select(diann_id, plink_id) %>%
                        na.omit %>%
                        unique
print(paste('Number of unique DIA-NN IDs:', length(unique(diannid_to_plinkid$diann_id))))
print(paste('Number of unique genetic IDs:', length(unique(diannid_to_plinkid$plink_id))))
# -

# ## Step 3 Map the protein IDs to plink IDs and set the file to contain columns: FID, IID and names of protein groups

proteins.out <- proteins %>%
                    merge(diannid_to_plinkid, all.x=T) %>%
                    select(FID=plink_id, IID=plink_id, everything()) %>%
                    subset(!is.na(IID))

# ## Step 4 Subset to only European ancestry

ancestry %<>% read.csv(sep='\t') %>% subset(Cluster=='EUR')

# Check how many out of our 1,551 individuals are EUR ancestry
proteins.out %>% merge(ancestry) %>% pull(Cluster) %>% table

proteins.out %<>% merge(ancestry) %>%
                    select(FID, IID, everything()) %>%
                    select(-diann_id, -Cluster)

# ## Step 5 Write out the proteins dataset and mapping file

write.table(proteins.out, file=parser$output, row.names=F, quote=F, sep='\t')
write.table(mapping_file, file=parser$mapping_file, row.names=F, quote=F, sep='\t')


