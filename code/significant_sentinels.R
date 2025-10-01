library(readr)
library(dplyr)

sentinels <- read_tsv("all_sentinels_merged.txt")

sig_sentinels <- sentinels %>%
  filter(P > 7.30103)   # since your sentinel file has "P" = log10(p)

# Save filtered list
write_tsv(sig_sentinels, "significant_sentinels.txt")
