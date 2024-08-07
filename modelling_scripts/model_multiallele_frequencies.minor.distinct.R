rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)
require(relaimpo)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

major_df <- fread("results/allele_frequency_out/all_sites.major.csv")
df <- fread("results/allele_frequency_out/all_sites.minor.csv") %>%
  left_join(hookup %>% dplyr::select(mutation_name, ref_AA, var_AA)) %>%
  filter()

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

# Check if all monthly frequencies have been calculated
all(df$mutation_name %in% parsed_agg$mutation_name)

df %>%
  filter(!(mutation_name %in% parsed_agg$mutation_name))

# Get blosum62 scores
data(BLOSUM62)
changes <- df %>%
  distinct(ref_AA, var_AA)

change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 42
  row <- changes[i, ]
  
  row %>%
    mutate(blosum62_score = BLOSUM62[row$ref_AA, ifelse(row$var_AA == "stop", "*", row$var_AA)])
}

change_df <- bind_rows(change_morsels)

merged <- df %>%
  left_join(parsed_agg) %>%
  left_join(change_df) %>%
  filter(!grepl("stop", mutation_name)) %>%
  # Set mutations that were not found in monthly counts to zero
  mutate(global_n = replace_na(global_n, 0),
         max_prop = replace_na(max_prop, 0)) %>%
  mutate(is_fixed = max_prop > 0.5) %>%
  filter(global_n != 0) %>%
  mutate(norm_n = (n - mean(n)) / sd(n)) %>%
  mutate(norm_sd_freq = replace_na(norm_sd_freq, 0))

# fwrite(plot_df, "results/allele_frequency_out/freq_stats.csv")

table(merged$is_fixed)

merged %>%
  distinct(mutation_name) %>%
  nrow()

## Explore variables ##
# Global frequency
merged %>%
  ggplot(aes(x = max_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/max_freq_versus_global_n.significant.minor.png", width = 4, height = 4)

merged %>%
  ggplot(aes(x = median_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Median intrahost freq.", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/median_freq_versus_global_n.significant.minor.png", width = 4, height = 4)

merged %>%
  ggplot(aes(x = norm_n, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "No. biosamples detected", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/n_versus_global_n.siginficant.minor.png", width = 4, height = 4)

merged %>%
  ggplot(aes(x = factor(blosum62_score), y = log10(global_n), fill = factor(blosum62_score))) +
  geom_boxplot() +
  # scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "BLOSUM62 score", y = "Log10(no. GISAID sequences)") +
  scale_fill_viridis_d()

ggsave("results/allele_frequency_out/blosum62_versus_global_n.significant.minor.png", width = 4, height = 4)

merged %>%
  ggplot(aes(x = norm_sd_freq, y = log10(global_n))) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", color = "black", fill = "orange") +
  theme_bw() +
  labs(x = "Coefficient of variance (allele freq.)", y = "Log10(no. GISAID sequences)")

ggsave("results/allele_frequency_out/CoV_versus_global_n.significant.minor.png", width = 4, height = 4)

# Model global counts
linreg <- lm(log(global_n) ~ norm_n + norm_sd_freq + max_freq + blosum62_score, 
             data = merged)

summary(linreg)
res <- anova(linreg)
res

prop_var <- calc.relimp(linreg)$lmg * 100

as.data.frame(coefficients(summary(linreg))) %>%
  rownames_to_column("predictor") %>%
  left_join(tibble(predictor = names(prop_var),
                   perc_var_expl = prop_var) %>%
              arrange(desc(perc_var_expl)))

# GLM
qpois <- glm(global_n ~ norm_n + norm_sd_freq + max_freq + blosum62_score,
             family = quasipoisson,
             data = merged)
coefficients(summary(qpois))

aod <- data.frame(anova(qpois))
null_deviance <- aod["NULL", "Resid..Dev"]

aod %>%
  dplyr::select(Deviance) %>%
  mutate(dev_explained = Deviance / null_deviance * 100) 

(null_deviance - qpois$deviance) / null_deviance * 100

# Monthly frequency
merged %>%
  ggplot(aes(x = is_fixed, y = max_freq)) +
  geom_boxplot() +
  geom_pwc()

merged %>%
  ggplot(aes(x = is_fixed, y = blosum62_score)) +
  geom_boxplot() +
  geom_pwc()

merged %>%
  ggplot(aes(x = is_fixed, y = norm_n)) +
  geom_boxplot() +
  geom_pwc()

merged %>%
  group_by(is_fixed) %>%
  summarise(n = n())
ggplot(aes(x = is_fixed, y = norm_sd_freq)) +
  geom_boxplot() +
  geom_pwc()

# 
merged %>%
  ggplot(aes(x = norm_n, y = median_prop)) +
  geom_point()


merged %>%
  mutate(is_significant = max_prop > 0.05) %>%
  group_by(is_significant) %>%
  summarise(n = n())

bind_rows(morsels) %>%
  ggplot(aes(freq_threshold, prop)) +
  geom_line()
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
# # 
# logreg <- glm(is_fixed ~ norm_sd_freq + norm_n + median_freq + max_freq + blosum62_score,
#               data = merged,
#               family = "binomial")
# # 
# anova(logreg)
# summary(logreg)
# # 
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
