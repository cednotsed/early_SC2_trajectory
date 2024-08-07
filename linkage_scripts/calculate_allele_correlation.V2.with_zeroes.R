rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

freq_df <- fread("results/allele_frequency_out/all_sites.raw.csv") %>%
  filter(freq > 0.02)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

mat <- freq_df %>%
  group_by(id, mutation_name) %>%
  summarise(sum_freq = sum(freq)) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = mutation_name, values_from = sum_freq) %>%
  column_to_rownames("id")

mat[is.na(mat)] <- 0

# Get mutation counts
col_sums <- colSums(mat > 0)
count_df <- tibble(mutation = names(col_sums), 
                   n_biosamples = col_sums)
to_keep <- names(col_sums)[col_sums > 5]

mat_filt <- mat[, to_keep]

pairs <- combn(to_keep, 2)

print(ncol(pairs))

morsels <- foreach(i = seq(ncol(pairs))) %do% {
  print(i)
  pair <- c("N_R203K", "N_G204R")
  pair <- pairs[, i]
  # pair <- c("NS3_G251V", "N_Q239H")
  mat_temp <- mat[, pair]
  mat_bool <- mat_temp > 0
  samples_to_keep <- rowSums(mat_bool) > 1
  samples_to_keep2 <- rowSums(mat_bool) == 2
  
  if(sum(samples_to_keep) >= 5) {
    # Remove double zeroes
    mat_temp <- mat_temp[samples_to_keep, ]
    corr <- cor(mat_temp[, 1], mat_temp[, 2])
    mut1_max <- max(mat_temp[, 1])
    mut2_max <- max(mat_temp[, 2])
    
    # Remove single zeroes
    mat_temp <- mat_temp[samples_to_keep, ]
    corr_no_zeroes <- cor(mat_temp[, 1], mat_temp[, 2])
    
    tibble(mutation1 = pair[1], mutation2 = pair[2], 
           corr = corr, 
           corr_no_zeroes = corr_no_zeroes,
           n_comparisons = sum(samples_to_keep),
           mut1_max = mut1_max,
           mut2_max = mut2_max)
  } else {
    return(NULL)
  }
}

merged <- bind_rows(morsels) %>%
  arrange(desc(corr)) %>%
  mutate(rsquared = corr * corr) %>%
  left_join(count_df %>% dplyr::rename(mutation1 = mutation, n_biosamples1 = n_biosamples)) %>%
  left_join(count_df %>% dplyr::rename(mutation2 = mutation, n_biosamples2 = n_biosamples)) %>%
  filter(!is.na(corr))

merged_filt <- merged %>%
  filter(rsquared > 0.25)

fwrite(merged, "results/linkage_out/linkage.n5.with_zeroes.csv")
fwrite(merged_filt, "results/linkage_out/linkage.n5.with_zeroes.rsquared25.csv")

