rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)

snp_df <- fread("results/pipeline_out/snp_counts.csv") %>%
  filter(snps <= 10)

fna <- readDNAStringSet("data/alignments/reassembled/reassembled.aln")
fna <- fna[names(fna) %in% snp_df$genome_name]

matrix_to_dnaset <- function(mat) {
  temp_mat <- apply(mat, 1, paste0, collapse = "")
  dnaset <- DNAStringSet(temp_mat, use.names = T)
  return(dnaset)
}

mat <- as.matrix(fna)

# Remove gappy sequences
prop_gaps <- apply(mat, 1,
                   function(x) {sum(x %in% c("-", "N")) / ncol(mat)})

to_remove <- prop_gaps > 0.10
n_removed <- length(prop_gaps[to_remove])

print(str_glue("Removed {n_removed} sequences"))

mat <- mat[!to_remove, ]

# Mask gappy sites
prop_site_gaps <- apply(mat, 2,
                        function(x) {sum(x %in% c("-", "N")) / nrow(mat)})

site_to_mask <- prop_site_gaps > 0.10

mat[, site_to_mask] <- "N"

masked_aln <- matrix_to_dnaset(mat)

writeXStringSet(masked_aln, "data/alignments/reassembled/reassembled.masked.aln")
