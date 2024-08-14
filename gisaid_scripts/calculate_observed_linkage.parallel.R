rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

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
print(ncol(pairs))
pair_index <- seq(ncol(pairs))

index_list <- split(pair_index, ceiling(seq_along(pair_index) / 10000))

cl <- makeCluster(12)
registerDoParallel(cl)

foreach(idx = seq(length(index_list)), 
        .packages = c("tidyverse", "data.table", "foreach")) %dopar% {
  temp_idx_list <- index_list[[idx]]
  chunk <- pairs[, temp_idx_list]
  
  morsels <- foreach(i = seq(ncol(chunk))) %do% {
    # print(i)
    pair <- chunk[, i]
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
    fwrite(str_glue("results/linkage_out/observed_linkage.temp/observed_linkage.gt1000.{idx}.csv"))
  
  return(NULL)
}
