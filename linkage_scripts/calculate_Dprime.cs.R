rm(list = ls())
setwd("/SAN/ugi/HAP_VAP/early_SC2_trajectory")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

meta <- fread("data/metadata/all_sra_metadata.csv")
timeframes <- meta %>%
  group_by(dataset, alias) %>%
  summarise(start = min(collection_month),
            end = max(collection_month)) %>%
  mutate(start = start %m-% months(1),
         end = end %m+% months(1))

args <- commandArgs(trailingOnly=TRUE)
d <- args[1]

timeframe_filt <- timeframes %>% 
  filter(alias == d)

monthly_df <- fread("results/allele_frequency_out/observed_out/all_monthly_frequencies.csv")

parsed_agg <- fread(str_glue("results/allele_frequency_out/observed_out/before_{d}.aggregate.csv"))

# Get total genome count
monthly_filt <- monthly_df %>%
  filter(collection_month <= timeframe_filt$start) %>%
  distinct(collection_month, n_total)

total_n <- sum(monthly_filt$n_total)

parsed_filt <- parsed_agg %>%
  filter(global_n > 1000)

savs <- parsed_filt$mutation_name

mat <- fread("data/metadata/gisaid/presence_absence_matrix.csv") %>%
  as_tibble()

colnames(mat) <- gsub("\\/", "", colnames(mat))

# Retain sequences within timeframe of dataset
mat_filt <- mat[mat$collection_month < timeframe_filt$start, colnames(mat) %in% savs]

# Total number of seqs. considered
total_n <- nrow(mat_filt)

dim(mat_filt)

pairs <- combn(colnames(mat_filt), 2)
print(ncol(pairs))

# Get SAV frequencies
freq_morsels <- foreach(mut = colnames(mat_filt)) %do% {
  # mut = colnames(mat_filt)[1]
  print(mut)
  p <- sum(mat_filt[, mut]) / total_n
  not_p <- sum(!mat_filt[, mut]) / total_n
  tibble(mut = mut, p = p, not_p = not_p)
}

sav_freq_df <- bind_rows(freq_morsels)

# Get pAB
link_morsels <- foreach(i = seq(ncol(pairs))) %do% {
  print(i)
  mut1 <- pairs[1, i]
  mut2 <- pairs[2, i]
  # print(paste(mut1, mut2))
  
  row_counts <- rowSums(mat_filt[, c(mut1, mut2)])
  pAB <- sum(row_counts == 2) / total_n
  tibble(mutation1 = mut1, mutation2 = mut2, pAB = pAB)
}

link_df <- bind_rows(link_morsels) %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
         mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
  left_join(sav_freq_df %>% dplyr::rename(mut1 = mut, pA = p, not_pA = not_p)) %>%
  left_join(sav_freq_df %>% dplyr::rename(mut2 = mut, pB = p, not_pB = not_p)) %>%
  mutate(D = pAB - pA * pB) %>%
  rowwise() %>%
  mutate(Dmax = ifelse(D > 0, 
                       min(c(pA * (1 - pB), (1 - pA) * pB)),
                       min(c(pA * pB, (1 - pA) * (1 - pB))))) %>%
  mutate(Dprime = D / Dmax) %>%
  mutate(r2 = (D * D) / (pA * (1 - pA) * pB * (1 - pB)))

# Remove same codon mutations
hookup_df <- tibble(mutation_name = unique(c(link_df$mut1, link_df$mut2))) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
  select(mutation_name, codon_name)

link_filt <- link_df %>%
  left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
  left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
  filter(codon1 != codon2) %>%
  select(mut1, mut2, pAB, pA, pB, D, Dmax, Dprime, r2)

link_filt %>%
  fwrite(str_glue("results/linkage_out/observed_linkage/observed_Dprime.before_{d}.gt1000.csv"))

# hist(link_filt$r2)
# 
# cor.test(link_filt$r2, link_filt$Dprime, method = "spearman")
# 
# test <- fread("results/linkage_out/observed_linkage.before_alpha.gt1000.csv")
# 
# test %>%
#   select(mut1 = mutation1, mut2 = mutation2, corr) %>%
#   right_join(link_df %>% select(mut1, mut2, Dprime, r2)) %>%
#   ggplot(aes(abs(Dprime), r2)) +
#   geom_bin2d() +
#   geom_smooth(method = "lm")
# 
# test %>%
#   select(mut1 = mutation1, mut2 = mutation2, corr) %>%
#   right_join(link_df %>% select(mut1, mut2, Dprime, r2)) %>%
#   ggplot(aes(abs(corr), r2)) +
#   geom_bin2d() +
#   geom_smooth(method = "lm")
# 
# 
# link_df %>% filter(Dprime > 0)

