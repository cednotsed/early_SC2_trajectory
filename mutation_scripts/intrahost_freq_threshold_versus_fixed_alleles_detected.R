rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

df <- fread("results/allele_frequency_out/all_sites.significant.csv") %>%
  left_join(hookup %>% dplyr::select(mutation_name, ref_AA, var_AA))

# Model allele frequencies versus observed 'fitness'
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.distinct.csv")
monthly_agg <- monthly_df %>%
  filter(!(collection_month %in% c("2019-12", "2020-01", "2020-02", "2020-03"))) %>%
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

freq_df <- fread("results/allele_frequency_out/all_sites.csv")

morsels <- foreach(t = seq(0, 1, 0.01)) %do% {
  freq_filt <- freq_df %>% filter(max_freq > t)  
  
  total_muts <- n_distinct(freq_filt$mutation_name)
  
  all <- monthly_df %>%
    group_by(mutation_name) %>%
    summarise(is_fixed = any(prop > 0.9)) %>%
    filter(is_fixed)
  
    n_fixed <- sum(all$mutation_name %in% freq_filt$mutation_name)
  
  tibble(freq_threshold = t, 
         n_fixed = n_fixed,
         n_total = total_muts,
         prop = n_fixed / nrow(all))
}

bind_rows(morsels) %>% 
  ggplot(aes(x = n_fixed / n_total, y = prop, color = freq_threshold)) +
  geom_point() +
  theme_bw() + 
  geom_vline(xintercept = 0.01, lty = "dashed") +
  labs(x = "Prop. of intrahost alleles", y = "Prop. of fixed alleles detected",
       color = "Freq. threshold")

ggsave("results/mutation_out/threshold_versus_fixed_prop.png", 
       dpi = 600, width = 4, height = 4)
