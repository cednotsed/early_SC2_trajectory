rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
pair_df <- fread("results/linkage_out/candidate_linkage.csv")
# 
# parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")
# 
# parsed_filt <- parsed_agg %>%
#   filter(global_n > 1000)

# # Master df
# cor_df <- fread("results/linkage_out/linkage.n5.with_zeroes.csv") %>%
#   mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
#          mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
#   select(n_comparisons, n_biosamples1, n_biosamples2, mut1, mut2, intra_corr = corr, corr_no_zeroes)
# 
freq_df <- fread("results/allele_frequency_out/all_sites.raw.csv") %>%
  filter(freq > 0.02)
# 
# pair_df <- cor_df %>%
#   filter(intra_corr * intra_corr > 0.1) %>%
#   filter(mut1 %in% parsed_filt$mutation_name, mut2 %in% parsed_filt$mutation_name) %>%
#   distinct(mut1, mut2)
# 
# hookup_df <- tibble(mutation_name = unique(c(pair_df$mut1, pair_df$mut2))) %>%
#   separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
#   mutate(codon_number = parse_number(mut)) %>%
#   mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
#   select(mutation_name, codon_name)
# 
# pair_filt <- pair_df %>%
#   left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
#   left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
#   filter(codon1 != codon2)

pdf("results/linkage_out/correlation_plots.with_zeroes.filt.pdf", width = 4, height = 4)
for(i in seq(nrow(pair_df))) {
  # i = 3
  row <- pair_df[i, ]
  pair <- str_glue("{row$mut1}-{row$mut2}")
  
  temp_df <- freq_df %>%
    filter(mutation_name %in% c(row$mut1, row$mut2)) %>%
    pivot_wider(id_cols = id, names_from = mutation_name, values_from = freq) %>% 
    mutate_if(is.numeric, ~replace_na(., 0))
  
  corr <- cor.test(temp_df[[row$mut1]], temp_df[[row$mut2]])
  r <- signif(corr$estimate, 2)
  p <- signif(corr$p.value, 2)
  
  temp_plot <- temp_df %>%
    ggplot(aes(x = get(row$mut1), y = get(row$mut2))) +
    geom_point(fill = "darkviolet", color = "black", alpha = 0.5,
               size = 3,
               pch = 21) +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    theme_bw() +
    theme(text = element_text(family = "sans"),
          axis.title = element_text(face = "bold")) +
    labs(x = row$mut1, y = row$mut2) +
    annotate("text", x = 0, y = 1, label = str_glue("r = {r}, p = {p}"),
             hjust = 0)
  print(temp_plot) 
}
dev.off()
