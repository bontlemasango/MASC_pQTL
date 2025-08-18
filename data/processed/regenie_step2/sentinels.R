sink(stderr())

args <- commandArgs(TRUE)
argc <- length(args)

WIDTH <- 1000000  # default 1Mb

if (argc < 1) {
  cat("Usage: sentinels.R HITS [ WIDTH (default", WIDTH, ")]\n")
  q()
}

hits <- args[1]
if (argc > 1) {
  WIDTH <- as.integer(args[2])
}

suppressMessages(library(data.table))

# Function to select sentinel SNPs
selectSentinels <- function(data.dt) {
  regions.list <- list()
  while (nrow(data.dt) > 0) {
    top.dt <- data.dt[which.max(P)]
    regions.list[[top.dt$rsid]] <- top.dt
    data.dt <- data.dt[
      chromosome != top.dt$chromosome |
      position <= top.dt$position - WIDTH |
      position >= top.dt$position + WIDTH
    ]
  }
  rbindlist(regions.list)
}

results.dt <- fread(cmd = paste("cat", hits))
setnames(results.dt, "ID", "rsid")
setnames(results.dt, "CHROM", "chromosome")
setnames(results.dt, "GENPOS", "position")
setnames(results.dt, "LOG10P", "P")

sentinels.dt <- selectSentinels(results.dt)
setnames(sentinels.dt, c("chromosome", "position"), c("chrom", "pos"))
setkeyv(sentinels.dt, c("chrom", "pos"))

sink()
# Print number of SNPs selected
cat("Number of SNPs selected:", nrow(sentinels.dt), "\n")

# Write output to the file provided
write.table(sentinels.dt, file = args[3], row.names = FALSE, quote = FALSE, sep = "\t")



