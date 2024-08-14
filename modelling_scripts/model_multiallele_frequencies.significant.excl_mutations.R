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

df <- fread("results/allele_frequency_out/all_sites.significant.csv") %>%
  left_join(hookup %>% distinct(mutation_name, ref_AA, var_AA))

# Model allele frequencies versus observed 'fitness'
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
monthly_agg <- monthly_df %>%
  group_by(mutation_name, collection_month) %>%
  summarise(monthly_n = sum(n_present),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

parsed_agg <- monthly_agg %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            max_prop = max(monthly_prop)) %>%
  filter(mutation_name != "")

# Check if all monthly frequencies have been calculated
all(df$mutation_name %in% parsed_agg$mutation_name)

df %>%
  filter(!(mutation_name %in% parsed_agg$mutation_name))

# Get blosum62 scores
aa_df <- fread("data/metadata/aa_properties.csv")

data(BLOSUM62)
changes <- df %>%
  distinct(ref_AA, var_AA)

change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 42
  row <- changes[i, ]
  
  row %>%
    left_join(aa_df %>% dplyr::select(ref_AA = amino_acid, ref_mw = mw, ref_charge = charge, ref_hydropathy = hydropathy)) %>%
    left_join(aa_df %>% dplyr::select(var_AA = amino_acid, var_mw = mw, var_charge = charge, var_hydropathy = hydropathy)) %>%
    mutate(delta_mw = var_mw - ref_mw,
           delta_charge = var_charge - ref_charge,
           delta_hydropathy = var_hydropathy - ref_hydropathy) %>%
    mutate(abs_mw = abs(delta_mw),
           abs_charge = abs(delta_charge),
           abs_hydropathy = abs(delta_hydropathy)) %>%
    mutate(blosum62_score = BLOSUM62[row$ref_AA, ifelse(row$var_AA == "stop", "*", row$var_AA)])
}

change_df <- bind_rows(change_morsels)

monthly_df %>%
  filter(mutation_name %in% c("Spike_D614G", "NSP12_P323L", "N_G204R", "N_R203K")) %>%
  filter(collection_month == "2020-03")

merged <- df %>%
  left_join(parsed_agg) %>%
  left_join(change_df) %>%
  left_join(hookup %>% distinct(mutation_name, region, protein_name)) %>%
  filter(!grepl("stop", mutation_name)) %>%
  # Set mutations that were not found in monthly counts to zero
  mutate(global_n = replace_na(global_n, 0),
         max_prop = replace_na(max_prop, 0)) %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  mutate(norm_sd_freq = replace_na(norm_sd_freq, 0)) %>%
  filter(!grepl("\\*", mutation_name)) %>%
  filter(!(mutation_name %in% c("Spike_D614G", "NSP12_P323L")))

linreg <- lm(log10(global_n + 1) ~ n + max_freq + blosum62_score + abs_mw + abs_charge + abs_hydropathy, 
             data = merged)

# linreg <- lm(log10(global_n + 1) ~ n + norm_sd_freq + max_freq, 
#              data = merged)
summary(linreg)
res <- anova(linreg)

hist(linreg$residuals)
# plot(linreg)

# shapiro.test(linreg$residuals)
prop_var <- calc.relimp(linreg)$lmg * 100
as.data.frame(coefficients(summary(linreg))) %>%
  rownames_to_column("predictor") %>%
  left_join(tibble(predictor = names(prop_var),
                   perc_var_expl = prop_var) %>%
              arrange(desc(perc_var_expl))) %>%
  mutate_if(is.numeric, function(x) {signif(x, 2)}) %>%
  arrange(desc(perc_var_expl))

rlinreg <- MASS::rlm(log(global_n + 1) ~ n + norm_sd_freq + max_freq + blosum62_score + abs_mw + abs_charge + abs_hydropathy, 
                     data = merged)

summary(rlinreg)
res <- anova(rlinreg)

sfsmisc::f.robftest(rlinreg, var = "n")
sfsmisc::f.robftest(rlinreg, var = "norm_sd_freq")
sfsmisc::f.robftest(rlinreg, var = "max_freq")
sfsmisc::f.robftest(rlinreg, var = "blosum62_score")
sfsmisc::f.robftest(rlinreg, var = "abs_mw")
sfsmisc::f.robftest(rlinreg, var = "abs_charge")
sfsmisc::f.robftest(rlinreg, var = "abs_hydropathy")

# GLM
qpois <- glm(global_n ~ n + norm_sd_freq + max_freq + blosum62_score + abs_mw + abs_charge + abs_hydropathy,
             family = quasipoisson,
             data = merged)
coefficients(summary(qpois))

aod <- data.frame(anova(qpois))
null_deviance <- aod["NULL", "Resid..Dev"]
summary(qpois)
aod %>%
  dplyr::select(Deviance) %>%
  mutate(dev_explained = Deviance / null_deviance * 100) 

(null_deviance - qpois$deviance) / null_deviance * 100

# genadd <- gam(global_n ~ s(n) + s(n_filt) + n_alleles + s(median_freq) + s(max_freq), 
#               data=plot_df)
# 
# res <- summary(genadd)
# obs_dev_explained <- res$dev.expl
# obs_dev_explained

# draw(genadd, residuals = T)


# # Monthly frequency
merged %>%
  ggplot(aes(x = is_fixed, y = max_freq)) +
  geom_boxplot() +
  geom_pwc()

merged %>%
  ggplot(aes(x = norm_n, y = max_prop)) +
  geom_point()


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
# logreg <- glm(is_fixed ~ norm_sd_freq + norm_n + blosum62_score + median_freq + max_freq,
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
