rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  distinct(mutation_name, protein_name, codon_number)

obs_df <- fread("results/linkage_out/observed_linkage.csv")

mut_df <- tibble(mutation_name = unique(c(obs_df$mutation1, obs_df$mutation2))) %>%
  separate(mutation_name, c("protein_name", "mut"), remove = F) %>%
  mutate(codon_number = extract_numeric(mut)) %>%
  select(mutation_name, protein_name, codon_number)

index_df <- hookup %>%
  distinct(protein_name, codon_number) %>%
  filter(codon_number != -1)

index_df <- mut_df %>%
  left_join(index_df %>% mutate(protein_index = seq(nrow(index_df))))

merged <- obs_df %>%
  left_join(index_df %>% dplyr::select(mutation1 = mutation_name, index1 = protein_index)) %>%
  left_join(index_df %>% dplyr::select(mutation2 = mutation_name, index2 = protein_index)) %>%
  mutate(distance = abs(index1 - index2))

merged %>%
  filter(corr > 0) %>%
  # filter(corr > 0.5) %>%
  filter(corr > 0) %>%
  ggplot(aes(x = distance, y = corr)) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_smooth() +
  theme_bw()
  
test <- lm(corr ~ 1/distance,
   data = merged)


