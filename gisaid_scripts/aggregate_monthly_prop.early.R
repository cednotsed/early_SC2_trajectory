rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv") %>%
  filter(collection_month %in% c("2019-12", "2020-01", "2020-02",
                                 "2020-03", "2020-04", "2020-05",
                                 "2020-06", "2020-07", "2020-08",
                                 "2020-09", "2020-10", "2020-11"))

possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q")

monthly_agg <- monthly_df %>%
  group_by(mutation_name, collection_month) %>%
  summarise(monthly_n = sum(n_present),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

parsed_agg <- monthly_agg %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            median_prop = median(monthly_prop),
            max_prop = max(monthly_prop)) %>%
  filter(mutation_name != "") %>%
  filter(!grepl("del|ins|stop", mutation_name)) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(mut = gsub("[0-9]", "", mut)) %>% 
  filter(mut != "") %>% 
  # mutate(mut = ifelse(grepl("stop", mut),
  #                     gsub("stop", "*", mut),
  #                     mut)) %>%
  mutate(ref_AA = substr(mut, 1, 1),
         var_AA = substr(mut, 2, 2)) %>%
  filter(ref_AA %in% possible_variants & var_AA %in% possible_variants)

parsed_agg %>%
  fwrite("results/mutation_out/monthly_freq_aggregate.early.csv")

