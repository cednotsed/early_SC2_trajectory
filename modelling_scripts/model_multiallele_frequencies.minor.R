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

df <- fread("results/allele_frequency_out/all_sites.minor.csv") %>%
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

merged <- df %>%
  left_join(parsed_agg) %>%
  left_join(change_df) %>%
  filter(!grepl("stop", mutation_name)) %>%
  # Set mutations that were not found in monthly counts to zero
  mutate(global_n = replace_na(global_n, 0),
         max_prop = replace_na(max_prop, 0)) %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  mutate(norm_sd_freq = replace_na(norm_sd_freq, 0)) %>%
  filter(!grepl("\\*", mutation_name))

# fwrite(plot_df, "results/allele_frequency_out/freq_stats.csv")

table(merged$is_fixed)

merged %>%
  distinct(mutation_name) %>%
  nrow()

## Explore variables ##
# Global frequency
cor_test <- cor.test(merged$max_freq, merged$global_n, method = "spearman")
rho <- signif(cor_test$estimate, 2)

p1 <- merged %>%
  ggplot(aes(x = max_freq, y = log10(global_n))) +
  # geom_point(color = "coral") +
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_smooth(method = "lm", color = "black", fill = "grey") +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Max intrahost freq.", y = "Log10(no. GISAID sequences)",
       fill = "No. mutations") +
  annotate("text", x = 0.3, y = 7, label = str_glue("rho={rho}, p<0.0001"))

ggsave("results/allele_frequency_out/max_freq_versus_global_n.minor.pdf", 
       plot = p1, width = 4, height = 4)

# cor_test <- cor.test(merged$median_freq, merged$global_n, method = "spearman")
# rho <- signif(cor_test$estimate, 2)
# 
# p1 <- merged %>%
#   ggplot(aes(x = median_freq, y = log10(global_n))) +
#   geom_bin2d() +
#   # geom_point(color = "steelblue") +
#   geom_smooth(method = "lm", color = "black", fill = "grey") +
#   theme_bw() +
#   theme(legend.position = "top",
#         text = element_text(family = "sans"),
#         axis.title = element_text(face = "bold"),
#         legend.title = element_text(face = "bold")) +
#   scale_fill_viridis_c() +
#   labs(x = "Median intrahost freq.", 
#        y = "Log10(no. GISAID sequences)", 
#        fill = "No. mutations") +
#   annotate("text", x = 0.3, y = 7, label = str_glue("rho={rho}, p<0.0001"))
# 
# ggsave("results/allele_frequency_out/median_freq_versus_global_n.significant.pdf", plot = p1, width = 4, height = 4)

cor_test <- cor.test(merged$n, merged$global_n, method = "spearman")
rho <- signif(cor_test$estimate, 2)

p2 <- merged %>%
  ggplot(aes(x = n, y = log10(global_n))) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_smooth(method = "lm", color = "black", fill = "grey") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "No. biosamples detected", 
       y = "Log10(no. GISAID sequences)", 
       fill = "No. mutations") +
  annotate("text", x = 50, y = 7, label = str_glue("rho={rho}, p<0.0001"))

ggsave("results/allele_frequency_out/n_versus_global_n.minor.pdf", 
       plot = p2, width = 4, height = 4)

cor_test <- cor.test(merged$norm_sd_freq, merged$global_n, method = "spearman")
rho <- signif(cor_test$estimate, 2)

p3 <- merged %>%
  ggplot(aes(x = norm_sd_freq, y = log10(global_n))) +
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_smooth(method = "lm", color = "black", fill = "grey") +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "Coefficient of variation", 
       y = "Log10(no. GISAID sequences)", 
       fill = "No. mutations") +
  annotate("text", x = 0, y = 7, label = str_glue("rho={rho}, p<0.0001"))

ggsave("results/allele_frequency_out/CoV_versus_global_n.minor.pdf", 
       plot = p3, width = 4, height = 4)

ggarrange(p1, p2, p3, nrow = 1)
ggsave("results/allele_frequency_out/sig_predictors.minor.pdf", 
       width = 12, height = 4)
# merged %>%
#   ggplot(aes(x = n_filt, y = log10(global_n))) +
#   geom_point(color = "steelblue") +
#   geom_smooth(color = "black", fill = "orange") +
#   theme_bw() +
#   labs(x = "No. biosamples detected (freq>0.1)", y = "Log10(no. GISAID sequences)")
# 
# ggsave("results/allele_frequency_out/n_filt_versus_global_n.png", width = 4, height = 4)
test1 <- cor.test(merged$median_freq, merged$global_n, method = "spearman")
test2 <- cor.test(merged$max_freq, merged$global_n, method = "spearman")
test3 <- cor.test(merged$n, merged$global_n, method = "spearman")
test4 <- cor.test(merged$norm_sd_freq, merged$global_n, method = "spearman")
test5 <- cor.test(merged$abs_charge, merged$global_n, method = "spearman")
test6 <- cor.test(merged$abs_mw, merged$global_n, method = "spearman")
test7 <- cor.test(merged$abs_hydropathy, merged$global_n, method = "spearman")
test8 <- cor.test(merged$blosum62_score, merged$global_n, method = "spearman")
tibble(var = c("median freq", "max freq", "n", "CV", "charge", "mw", "hydropathy", "blosum62"),
       rho = c(test1$estimate, test2$estimate, test3$estimate,
               test4$estimate, test5$estimate, test6$estimate,
               test7$estimate, test8$estimate),
       pval = c(test1$p.value, test2$p.value, test3$p.value,
                test4$p.value, test5$p.value, test6$p.value,
                test7$p.value, test8$p.value)) %>%
  mutate(rho = signif(rho, 2), 
         padj = signif(p.adjust(pval, method = "BH"), 2)) %>% View()

# merged %>%
#   ggplot(aes(x = factor(blosum62_score), y = log10(global_n), fill = factor(blosum62_score))) +
#   geom_boxplot() +
#   # scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "BLOSUM62 score", y = "Log10(no. GISAID sequences)")
# 
# ggsave("results/allele_frequency_out/blosum62_versus_global_n.significant.png", width = 4, height = 4)
# 
# merged %>%
#   ggplot(aes(x = abs_hydropathy, y = log10(global_n))) +
#   geom_point() +
#   geom_smooth(method = "lm")
#   # scale_fill_manual(values = c("olivedrab", "olivedrab3", "olivedrab1")) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "BLOSUM62 score", y = "Log10(no. GISAID sequences)")

# Model global counts
linreg <- lm(log10(global_n + 1) ~ n + max_freq + norm_sd_freq + blosum62_score + abs_mw + abs_charge + abs_hydropathy, 
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
