rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")
mat <- as.matrix(readDNAStringSet("data/alignments/reassembled/reassembled.aln"))

# Compute distances
ref_vec <- mat[1, ]

morsels <- foreach(i = seq(nrow(mat))) %do% {
  temp <- rbind(ref_vec, mat[i, ])
  not_Ns <- apply(temp, 2, function(x){all(x != "N")})
  temp_filt <- temp[, not_Ns]
  discordant <- apply(temp_filt, 2, function(x) {x[1] != x[2]})
  snp_count <- sum(discordant)
  tibble(genome_name = rownames(mat)[i], snps = snp_count)
}

merged <- bind_rows(morsels)

merged %>%
  fwrite("results/pipeline_out/snp_counts.csv")


merged %>%
  dplyr::rename(biosample = genome_name) %>%
  left_join(meta) %>% 
  filter(snps > 10) %>% View()
  ggplot(aes(x = factor(snps), y = as.Date(collection_date, "%Y-%m-%d"))) +
  geom_point()

merged %>%
  filter(snps > 10)
merged %>%
  arrange(desc(snps)) %>% View()
bind_rows(morsels) %>%
  ggplot(aes(x = snps)) +
  geom_histogram(bins = 100) +
  labs(x = "SNPs", y = "No. biosamples")

bind_rows(morsels) %>%
  filter(snps > 10)
  filter(genome_name == "SAMN13922059")
