rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)
require(igraph)

df <- fread("results/linkage_out/obs_cor.edgelist.csv")

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
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
  filter(mutation_name != "")

community = c(15)

muts <- df %>%
  distinct(mutation1, comm) %>%
  filter(comm %in% community)

monthly_agg %>%
  filter(mutation_name %in% muts$mutation1) %>%
  mutate(collection_month = as.Date(str_glue("{collection_month}-01"))) %>%
  ggplot(aes(x = collection_month, y = monthly_prop, color = mutation_name)) +
  geom_line() +
  geom_vline(xintercept = c(as.Date("2021-03-01"), as.Date("2021-06-01"),
                            as.Date("2022-01-01"), as.Date("2023-02-01"),
                            as.Date("2023-12-01")),
             lty = "dashed", color = "black")

