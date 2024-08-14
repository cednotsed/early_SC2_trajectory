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
  pair <- pairs[, i]
  # pair <- c("NS3_G251V", "N_Q239H")
  mat_temp <- mat[, pair]
  mat_bool <- mat_temp > 0
  samples_to_keep <- rowSums(mat_bool) == 2
  if(sum(samples_to_keep) > 5) {
    mat_temp <- mat_temp[samples_to_keep, ]
    corr <- cor(mat_temp[, 1], mat_temp[, 2])
    tibble(mutation1 = pair[1], mutation2 = pair[2], corr = corr, n_comparisons = sum(samples_to_keep))
  } else {
    return(NULL)
  }
}

merged <- bind_rows(morsels) %>%
  arrange(desc(corr)) %>%
  mutate(rsquared = corr * corr) %>%
  left_join(count_df %>% dplyr::rename(mutation1 = mutation, n_biosamples1 = n_biosamples)) %>%
  left_join(count_df %>% dplyr::rename(mutation2 = mutation, n_biosamples2 = n_biosamples))

merged_filt <- merged %>%
  filter(rsquared > 0.25)

fwrite(merged, "results/linkage_out/linkage.n5.csv")
fwrite(merged_filt, "results/linkage_out/linkage.n5.rsquared25.csv")

# gisaid_meta <- fread("data/metadata/gisaid/gisaid_metatadata.080724.filt.tsv")
# 
# gisaid_filt <- gisaid_meta %>%
#   filter(!(collection_month %in% c("2019-12", "2020-01", "2020-02", "2020-03")))
# 
# res_morsels <- foreach(i = seq(nrow(merged_filt))) %do% {
#   print(i)
#   row <- merged_filt[i, ]
#   mut1 <- row$mutation1
#   mut2 <- row$mutation2
#   
#   res <- gisaid_filt %>%
#     summarise(sum_linked = sum(grepl(mut1, aa_substitutions) & grepl(mut2, aa_substitutions)),
#               sum_unlinked = sum((grepl(mut1, aa_substitutions) & !grepl(mut2, aa_substitutions)) | 
#                                    (grepl(mut2, aa_substitutions) & !grepl(mut1, aa_substitutions))))
# 
#   temp <- row %>%
#     mutate(sum_linked = res$sum_linked,
#            sum_unlinked = res$sum_unlinked)
#   temp %>%
#     mutate(ratio = sum_linked / sum_unlinked) %>%
#     select(ratio)
#   return(temp)
# }

