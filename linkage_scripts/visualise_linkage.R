rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.distinct.csv")
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
  filter(global_n > 1)

# Master df
cor_df <- fread("results/linkage_out/linkage.n5.csv")

file_dir <- "results/linkage_out/monthly_linkage/"

file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

merged <- bind_rows(morsels) %>% 
  filter(mutation1 %in% parsed_filt$mutation_name,
         mutation2 %in% parsed_filt$mutation_name) %>%
  left_join(parsed_agg %>% dplyr::select(mutation1 = mutation_name, protein_name1 = protein_name, codon_number1 = codon_number)) %>%
  left_join(parsed_agg %>% dplyr::select(mutation2 = mutation_name, protein_name2 = protein_name, codon_number2 = codon_number)) %>%
  filter(!(protein_name1 == protein_name2 & codon_number1 == codon_number2)) %>%
  filter(!grepl("stop", mutation1),
         !grepl("stop", mutation2)) %>%
  group_by(mutation1, mutation2, corr, rsquared, n_biosamples1, n_biosamples2, n_comparisons) %>%
  summarise(n_linked = sum(n_linked),
            n_unlinked = sum(n_unlinked)) %>%
  mutate(ratio = ifelse(corr > 0, n_linked / n_unlinked, n_unlinked / n_linked)) 

obs_df <- fread("results/linkage_out/observed_linkage.csv") %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2)) %>%
  mutate(mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1))

merged %>% 
  ungroup() %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2)) %>%
  mutate(mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
  select(mut1, mut2, intra_corr = corr) %>%
  left_join(obs_df %>% select(mut1, mut2, obs_corr = corr)) %>% View()
  ggplot(aes(x = intra_corr, obs_corr)) +
  geom_point()

obs_df %>%
  filter(mutation1 == "E_L37F")

merged %>% View()
  filter(corr < 0) %>%
  arrange(desc(ratio)) %>% 
  View()

merged %>%
  filter(corr > 0) %>%
  arrange(ratio) %>% 
  ggplot(aes(x = corr, y = log10(ratio))) + 
  geom_point() +
  geom_smooth(method = "lm")
  
merged %>% View()
  
freq_df <- fread("results/allele_frequency_out/all_sites.raw.csv")

lr <- lm(log10(ratio + 0.00000000000001) ~ n_biosamples1 + n_biosamples2 + n_linked + n_unlinked + corr,
   data = merged %>% filter(corr > 0, 
                            ratio != Inf))

anova(lr)

merged %>%
  filter(corr > 0) %>%
  arrange(desc(corr))
summary(lr)
hist(lr$residuals)
mat <- freq_df %>%
  group_by(id, mutation_name) %>%
  summarise(sum_freq = sum(freq)) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = mutation_name, values_from = sum_freq) %>%
  column_to_rownames("id")

mat[is.na(mat)] <- 0
  
plot(mat[, "NSP12_K780R"], mat[, "NSP15_D91G"])

filter(corr > 0) %>% 
  arrange(desc(ratio)) %>% 
  filter(ratio != Inf) %>%
  filter(n_linked > 10) %>%
  ggplot(aes(x = corr, y = ratio)) +
  geom_point()
