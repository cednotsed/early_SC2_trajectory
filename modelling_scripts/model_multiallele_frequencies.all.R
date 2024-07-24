rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

fna <- readDNAStringSet("data/alignments/reassembled/reassembled.masked.aln")
meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  select(pos = nucleotide_pos, ref_nuc, var_nuc, 
         region, protein_name, ref_AA, 
         codon_number, pos_in_codon, codon_from_gene_start, var_AA, 
         mutation_type, mutation_name)

non_nanopore <- meta %>%
  filter(platforms != "OXFORD_NANOPORE")

df <- fread("results/allele_frequency_out/merged_frequencies.csv")

# # Check all allele frequency entries are unique
# df %>%
#   filter(id %in% names(fna)) %>%
#   filter(id %in% non_nanopore$biosample) %>%
#   group_by(id, pos, allele) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

df_filt <- df %>%
  filter(id %in% names(fna)) %>%
  filter(id %in% non_nanopore$biosample) %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(var_nuc = gsub("p_", "", allele)) %>%
  left_join(hookup %>% select(ref_nuc, pos, var_nuc, protein_name, pos_in_codon, ref_AA, codon_number))

# Check if there are any double/triple mutations
codon_allele_count <- df_filt %>%
  filter(ref_nuc != var_nuc) %>%
  group_by(protein_name, codon_number) %>%
  group_by(id, protein_name, codon_number, ref_AA) %>%
  summarise(n_codon_pos = n_distinct(pos_in_codon),
            n_alleles = n()) %>%
  arrange(desc(n_codon_pos)) %>% 
  ungroup() 

# Deconvolute mutations
# Chuck codons that have more than one allele per position
triple_df <- codon_allele_count %>%
  filter(n_codon_pos == 3) %>%
  filter(n_alleles == n_codon_pos)

double_df <- codon_allele_count %>%
  filter(n_codon_pos == 2) %>%
  filter(n_alleles == n_codon_pos)

single_df <- codon_allele_count %>%
  filter(n_codon_pos == 1)

# Parse multi-allele codons
triple_morsels <- foreach(i = seq(nrow(triple_df))) %do% {
  row <- triple_df[i, ]
  temp <- row %>%
    left_join(df_filt) %>%
    filter(ref_nuc != var_nuc) %>%
    arrange(pos_in_codon)
  
  # Check that alleles are in normal order
  if(all(temp$pos_in_codon == c(1, 2, 3))) {
    var_AA <- paste0(temp$var_nuc, collapse = "")
    var_AA <- as.character(translate(DNAString(var_AA)))
    freq <- mean(temp$freq)
    
    final <- temp %>%
      distinct(id, protein_name, ref_AA, codon_number) %>%
      mutate(var_AA = var_AA, freq = freq)
    
    return(final)
  }
}

double_df %>%
  distinct(ref_nuc)
  summarise(n = n_distinct(id),
            n_filt = n_distinct(id[freq > 0.1]),
            median_freq = median(freq),
            max_freq = max(freq)) %>%
  filter(!grepl("\\*", ref_AA)) # Remove stop codon mutations

arrange(desc(n)) %>%
  ungroup() %>%
  
parsed2 <- parsed

parsed2 %>% distinct(mutation_name) %>% nrow()

fwrite(parsed2, "results/allele_frequency_out/parsed_multiallelic_sites.csv")

# Model allele frequencies versus observed 'fitness'
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
monthly_agg <- monthly_df %>%
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
allele_count_per_mutation <- parsed %>%
  group_by(mutation_name) %>%
  summarise(n_alleles = n_distinct(pos, allele))

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

# Final mutations considered
plot_df <- parsed_filt %>%
  left_join(change_df) %>%
  # filter(n > 1) %>%
  mutate(is_fixed = max_prop > 0.5)

fwrite(plot_df, "results/allele_frequency_out/freq_stats.csv")

table(plot_df$is_fixed)

plot_df %>%
  distinct(pos, allele) %>%
  nrow()

plot_df %>%
  distinct(mutation_name) %>%
  nrow()

cor.test(plot_df$max_freq, plot_df$global_n)
# cor.test(plot_)
## Explore variables ##
# Global frequency
plot_df %>%
  ggplot(aes(x = max_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/max_freq_versus_global_n.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = median_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Median intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/median_freq_versus_global_n.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = n, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "No. biosamples detected", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_versus_global_n.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = n_filt, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "No. biosamples detected (freq>0.1)", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_filt_versus_global_n.png", width = 4, height = 4)

plot_df %>%
  ggplot(aes(x = factor(n_alleles), y = log10(global_n), fill = factor(n_alleles))) +
  geom_boxplot() +
  geom_pwc() +
  scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "No. alleles observed", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_alleles_versus_global_n.png", width = 4, height = 4)

# Model global counts
preg <- glm(global_n ~ n + n_filt + n_alleles + median_freq + max_freq,
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
