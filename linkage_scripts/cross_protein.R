rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  distinct(mutation_name, protein_name, codon_number)

obs_df <- fread("results/linkage_out/observed_linkage.csv")

merged <- obs_df %>%
  separate(mutation1, c("protein_name1"), "_", remove = F) %>%
  separate(mutation2, c("protein_name2"), "_", remove = F) %>%
  mutate(cross_protein = protein_name1 != protein_name2) %>%
  mutate(is_linked = corr <= 0.9)

merged %>%
  group_by(cross_protein) %>%
  summarise(n = n())
plot_df <- merged %>%
  group_by(is_linked, cross_protein) %>%
  summarise(n = n()) %>%
  group_by(is_linked) %>%
  mutate(prop = n / sum(n))

mat <- plot_df %>%
  select(is_linked, cross_protein, n) %>%
  pivot_wider(id_cols = is_linked, names_from = cross_protein, values_from = n) %>%
  column_to_rownames("is_linked")

fisher.test(mat)
ggplot(aes(x = is_linked, y = prop,  fill = cross_protein)) +
  geom_bar(stat = "identity")
merged %>%
  
  ggplot(aes(x = distance, y = corr)) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_smooth() +
  theme_bw()
