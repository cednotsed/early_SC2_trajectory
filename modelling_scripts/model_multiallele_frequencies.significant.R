rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)
require(gratia)

fna <- readDNAStringSet("data/alignments/reassembled/reassembled.masked.aln")
meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>%
  filter(!grepl("element|loop", protein_name)) %>%
  select(pos = nucleotide_pos, ref_nuc, var_nuc, 
         region, protein_name, ref_AA, 
         codon_number, codon_from_gene_start, var_AA, 
         mutation_type)

non_nanopore <- meta %>%
  filter(platforms != "OXFORD_NANOPORE")

nanopore <- meta %>%
  filter(platforms == "OXFORD_NANOPORE")

df <- fread("results/allele_frequency_out/merged_frequencies.csv")

# Check all allele frequency entries are unique
# df %>%
#   filter(id %in% names(fna)) %>%
#   filter(id %in% non_nanopore$biosample) %>%
#   group_by(id, pos, allele) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

df_filt <- df %>%
  filter(id %in% names(fna)) %>%
  filter(id %in% non_nanopore$biosample) %>%
  mutate(pos = as.numeric(pos))

parsed <- df_filt %>%
  filter(id %in% names(fna)) %>%
  filter(id %in% non_nanopore$biosample) %>%
  filter(freq > 0.1) %>%
  group_by(pos, allele) %>%
  summarise(n = n_distinct(id),
            n_filt = sum(freq > 0.1),
            median_freq = median(freq),
            max_freq = max(freq)) %>%
  mutate(var_nuc = gsub("p_", "", allele)) %>%
  arrange(desc(n)) %>%
  left_join(hookup)

parsed2 <- parsed %>%
  ungroup() %>%
  filter(mutation_type == "NS") %>%
  mutate(protein_name = gsub(" \\(part 2\\)| \\(part 1\\)", "", protein_name)) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(var_AA = ifelse(var_AA == "*", "stop", var_AA)) %>%
  mutate(protein_name = case_when(protein_name == "S" ~ "Spike",
                                  protein_name %in% c("NS3a", "NSP3a", "NSP12a", "NS8a") ~ gsub("a", "", protein_name),
                                  T ~ protein_name)) %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}")) %>%
  filter(!grepl("\\*", ref_AA)) # Remove stop codon mutations 

parsed2 %>% distinct(mutation_name) %>% nrow()

# fwrite(parsed2, "results/allele_frequency_out/parsed_multiallelic_sites.csv")

# Model allele frequencies versus observed 'fitness'
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
monthly_agg <- monthly_df %>%
  filter(!(collection_month %in% c("2019-12", "2020-01", "2020-02", "2020-03"))) %>%
  group_by(mutation_name, collection_month) %>%
  summarise(monthly_n = sum(n_seqs),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

parsed_agg <- monthly_agg %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            max_prop = max(monthly_prop)) %>%
  filter(!grepl("\\*", mutation_name))

# Check if all monthly frequencies have been calculated
all(parsed2$mutation_name %in% parsed_agg$mutation_name)

# Generate final dataframe
parsed2 %>%
  group_by(mutation_name) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# Get no. alleles for each mutation
allele_count_per_mutation <- parsed2 %>%
  group_by(mutation_name) %>%
  summarise(n_alleles = n())

# Get most common allele for each mutation
morsels <- foreach(mut = unique(parsed2$mutation_name)) %do% {
  print(mut)
  parsed2 %>%
    filter(mutation_name == mut) %>%
    arrange(desc(n), desc(median_freq)) %>%
    head(1)
}

parsed_filt <- bind_rows(morsels) %>%
  left_join(parsed_agg) %>%
  left_join(allele_count_per_mutation)

# Check all mutations are unique now
parsed_filt %>% distinct(mutation_name) %>% nrow() == parsed_filt %>% nrow()

# Get blosum62 scores
data(BLOSUM62)
changes <- parsed2 %>%
  distinct(ref_AA, var_AA)

change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 42
  row <- changes[i, ]
  
  row %>%
    mutate(blosum62_score = BLOSUM62[row$ref_AA, ifelse(row$var_AA == "stop", "*", row$var_AA)])
}

change_df <- bind_rows(change_morsels)

plot_df <- parsed_filt %>%
  left_join(change_df) %>%
  # filter(n > 1) %>%
  mutate(is_fixed = max_prop > 0.5)

fwrite(plot_df, "results/allele_frequency_out/freq_stats.significant.csv")

plot_df %>%
  distinct(pos, allele) %>%
  nrow()

plot_df %>%
  distinct(mutation_name) %>%
  nrow()

table(plot_df$is_fixed)

## Explore variables ##
# Global frequency
plot_df %>%
  ggplot(aes(x = max_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/max_freq_versus_global_n.significant.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = median_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Median intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/median_freq_versus_global_n.significant.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = n, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "No. biosamples detected", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_versus_global_n.significant.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = n_filt, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "No. biosamples detected (freq>0.1)", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_filt_versus_global_n.significant.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = factor(n_alleles), y = log10(global_n), fill = factor(n_alleles))) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "No. alleles observed", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_alleles_versus_global_n.significant.png", width = 4, height = 4)

# Model global counts
preg <- glm(global_n ~ n + n_filt + n_alleles + median_freq + max_freq + blosum62_score,
            family = poisson(),
            data = plot_df)

coefficients(summary(preg))
# preg <- glm(global_n ~ max_freq,
#             family = poisson(),
#             data = plot_df)

aod <- data.frame(anova(preg))
null_deviance <- aod["NULL", "Resid..Dev"]

aod %>%
  select(Deviance) %>%
  mutate(dev_explained = Deviance / null_deviance * 100) 

aod %>%
  select(Deviance) %>%
  mutate(dev_explained = Deviance / null_deviance * 100) %>%
  summarise(sum(dev_explained, na.rm = T))

# genadd <- gam(global_n ~ s(n) + s(n_filt) + n_alleles + s(median_freq) + s(max_freq), 
#               data=plot_df)
# 
# res <- summary(genadd)
# obs_dev_explained <- res$dev.expl
# obs_dev_explained

# draw(genadd, residuals = T)


# # Monthly frequency
# plot_df %>%
#   ggplot(aes(x = is_fixed, y = max_freq)) +
#   geom_boxplot() 
# geom_point(color = "steelblue") +
#   geom_smooth(color = "black", fill = "orange") +
#   theme_bw() +
#   labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)")
# 
# ggsave("results/allele_frequency_out/max_freq_versus_max_prop.png", width = 4, height = 4)
# 
# 
# plot_df %>%
#   group_by(is_fixed, n_alleles) %>%
#   summarise(n = n()) %>%
#   group_by(is_fixed) %>%
#   mutate(prop = n / sum(n)) %>%
#   ggplot(aes(x = is_fixed, y = n_alleles, fill = prop)) +
#   geom_tile()
# 
# plot_df %>%
#   ggplot(aes(x = median_freq, y = max_prop)) +
#   geom_point(color = "steelblue") +
#   geom_smooth(color = "black", fill = "orange") +
#   theme_bw() +
#   labs(x = "Median intrahost freq.", y = "Log10(no. GISAID sequences)")
# 
# ggsave("results/allele_frequency_out/median_freq_versus_global_n.png", width = 4, height = 4)
# 
# plot_df %>%
#   ggplot(aes(x = factor(n_alleles), y = log10(global_n), fill = factor(n_alleles))) +
#   geom_boxplot() +
#   geom_pwc() +
#   scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "No. alleles observed", y = "Log10(no. GISAID sequences)")
# 
# ggsave("results/allele_frequency_out/n_alleles_versus_global_n.png", width = 4, height = 4)
# 
# 
# 
# logreg <- glm(is_fixed ~ n + median_freq + max_freq,
#               data = plot_df,
#               family = "binomial")
# 
# anova(logreg)
# summary(logreg)
# 
# exp(coefficients(summary(logreg))[,"Estimate"])
# 
# preg <- glm(global_n ~ n + n_filt + max_freq + median_freq + ,
#             family = poisson(),
#             data = plot_df)
# 
# 
# anova(preg)
# summary(preg)
# plot_df %>%
#   filter(n > 1) %>%
#   filter(global_n != 0) %>%
#   ggplot(aes(x = n, y = global_n)) +
#   geom_bin2d(bins = 100) +
#   scale_fill_viridis_c() +
#   theme_classic() +
#   labs(x = "No. biosamples with allele detected", y = "Max. monthly frequency", 
#        fill = "Density")
# 
# plot_df %>%
#   ggplot(aes(x = is_fixed, y = log10(max_freq))) +
#   geom_boxplot() +
#   geom_pwc()
# 
# parsed2 %>%
#   filter(median_freq > 0.99) %>% View()
