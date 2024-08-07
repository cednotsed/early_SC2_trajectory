rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
monthly_agg <- monthly_df %>%
  group_by(mutation_name, collection_month) %>%
  summarise(monthly_n = sum(n_present),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q",
                       "*")

parsed_agg <- monthly_agg %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            median_prop = median(monthly_prop),
            max_prop = max(monthly_prop)) %>%
  filter(mutation_name != "") %>%
  filter(!grepl("del|ins", mutation_name)) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(mut = gsub("[0-9]", "", mut)) %>% 
  filter(mut != "") %>% 
  mutate(mut = ifelse(grepl("stop", mut), 
                      gsub("stop", "*", mut), 
                      mut)) %>% 
  mutate(ref_AA = substr(mut, 1, 1),
         var_AA = substr(mut, 2, 2)) %>%
  # Remove non-canonical amino acids
  filter(ref_AA %in% possible_variants,
         var_AA %in% possible_variants)

parsed_filt <- parsed_agg %>% 
  filter(global_n > 1000)

mat <- monthly_df %>%
  filter(mutation_name %in% parsed_filt$mutation_name) %>%
  select(mutation_name, collection_month, prop) %>%
  pivot_wider(id_cols = mutation_name, names_from = collection_month, values_from = prop) %>%
  column_to_rownames("mutation_name")

mat[is.na(mat)] <- 0

pairs <- combn(unique(parsed_filt$mutation_name), 2)

# Prune comparisons
# tibble(mutation1 = t(pairs)[, 1], mutation2 = t(pairs)[, 2]) %>%
#   left_join(parsed_agg %>% dplyr::select(mutation1 = mutation_name, protein_name1 = protein_name, codon_number1 = codon_number)) %>%
#   left_join(parsed_agg %>% dplyr::select(mutation2 = mutation_name, protein_name2 = protein_name, codon_number2 = codon_number)) %>%
#   filter(!(protein_name1 == protein_name2 & codon_number1 == codon_number2)) %>% nrow()
print(ncol(pairs))

morsels <- foreach(i = seq(ncol(pairs))) %do% {
  print(i)
  pair <- pairs[, i]
  # pair <- c("NS3_G251V", "N_Q239H")
  mat_temp <- t(mat[pair, ])
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

bind_rows(morsels) %>% 
  arrange(desc(corr)) %>%
  fwrite("results/linkage_out/observed_linkage.gt1000.csv")
 