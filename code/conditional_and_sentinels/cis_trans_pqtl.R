#read in the hgnc, conditional_analysis, and Uniport IDs files
library(readr)
library(readxl)
library(dplyr)
library(writexl)
hgnc_ucsc_b38_genes <- read.delim("~/Downloads/pQTL/hgnc_ucsc_b38_genes.tsv")
conditional_analysis <- read.table("~/Downloads/pQTL/conditional_analysis.txt",header = TRUE, sep = "\t")
Proteins_Uniport_Ids <- read_excel("~/Downloads/pQTL/Proteins_Uniport_Ids.xlsx")

#rename the column in hgnc to make it the same as the one with uniport Ids
names(hgnc_ucsc_b38_genes)[names(hgnc_ucsc_b38_genes) == "symbol"] <- "Gene_Symbol"

#rename the uniport_ID column accordindly
Proteins_Uniport_Ids <- Proteins_Uniport_Ids %>% rename(uniprot_ids = UniProt_ID)
# Ensure uniprot_ids is character in both data frames
hgnc_ucsc_b38_genes$uniprot_ids <- as.character(hgnc_ucsc_b38_genes$uniprot_ids)
Proteins_Uniport_Ids$uniprot_ids <- as.character(Proteins_Uniport_Ids$uniprot_ids)

#merge to create a new DF:new_hgnc_2
new_hgnc_2 <- merge(hgnc_ucsc_b38_genes, Proteins_Uniport_Ids, by = "uniprot_ids")

#rename a column in new_hgnc
names(new_hgnc_2)[names(new_hgnc_2) == "name_in_MASC"] <- "PROTEIN_ID"

#change the proteins with - and _ to have a full stop (.)
new_hgnc_2 <- new_hgnc_2 %>%
  mutate(PROTEIN_ID = gsub("[-_]", ".", PROTEIN_ID))

#rename a column, Gene_Symbol.y to Gene_Symbol
new_hgnc_2 <- new_hgnc_2 %>%
  rename(Gene_Symbol = Gene_Symbol.y)

#remove certain columns in new_hgnc_2
new_hgnc_2  <- new_hgnc_2 %>% select(-c(4:11))
new_hgnc_2  <- new_hgnc_2 %>% select(-c(5:11))
new_hgnc_2  <- new_hgnc_2 %>% select(-c(7:45))
new_hgnc_2$gencc <- NULL

#rename the chr column
new_hgnc_2 <- new_hgnc_2 %>%
  rename(CHROM = chrom)

#remove the "chr" prefix from the CHROM column
new_hgnc_2$CHROM <- gsub("chr", "", new_hgnc_2$CHROM)
#remove the sufix
new_hgnc_2$CHROM <- gsub("_.*", "", new_hgnc_2$CHROM)
#replace "X" with "23" in the CHROM column 
new_hgnc_2$CHROM[new_hgnc_2$CHROM == "X"] <- "23"

names(conditional_analysis)[names(conditional_analysis) == "Protein"] <- "PROTEIN_ID"

# Merge on protein ID
merged_dataframe <- merge(conditional_analysis, new_hgnc_2, by = "PROTEIN_ID")

#To calculate cis vs trans pQTLs using the rule Brandon sent, based on the merged_dataframe, assuming:
#Chr → SNP chromosome (from conditional_analysis)
#bp → SNP position
#CHROM → Gene chromosome (from new_hgnc_2)
#chromStart → Gene position

# Ensure chromosomes are treated as character (to avoid factor issues)
merged_dataframe$Chr <- as.character(merged_dataframe$Chr)
merged_dataframe$CHROM <- as.character(merged_dataframe$CHROM)

# Create a new column to classify pQTLs as cis or trans
merged_dataframe$pQTL_type <- ifelse(
  merged_dataframe$Chr == merged_dataframe$CHROM &
    abs(merged_dataframe$bp - merged_dataframe$chromStart) <= 1e6,
  "cis",
  "trans"
)

#count how many are trans and cis
table(merged_dataframe$pQTL_type)

#To get the range of cis-pQTL sentinel distances (i.e., the distance between the sentinel SNP and the corresponding gene start for cis pQTLs)
# Filter for cis-pQTLs only
cis_sentinels <- merged_dataframe %>% filter(pQTL_type == "cis")

# Calculate absolute distance between SNP and gene start
cis_sentinels <- cis_sentinels %>%
  mutate(cis_distance = abs(bp - chromStart))

# Get the range of distances
range_cis <- range(cis_sentinels$cis_distance, na.rm = TRUE)

# Print the result
print(range_cis)

#checks how many unique proteins are there:
length(unique(merged_dataframe$PROTEIN_ID))

#save the df
write_xlsx(merged_dataframe, "~/Downloads/pQTL/MASC_pQTL_22_09_2025.xlsx")

########################

