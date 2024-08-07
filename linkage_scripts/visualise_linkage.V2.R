rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(ggrepel)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")

# Master df
cor_df <- fread("results/linkage_out/linkage.n5.with_zeroes.csv") %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
         mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
  select(n_comparisons, mut1, mut2, mut1_max, mut2_max, corr_no_zeroes, intra_corr = corr)

# Remove same codon mutations
hookup_df <- tibble(mutation_name = unique(c(cor_df$mut1, cor_df$mut2))) %>%
  separate(mutation_name, c("protein_name", "mut"), "_", remove = F) %>%
  mutate(codon_number = parse_number(mut)) %>%
  mutate(codon_name = str_glue("{protein_name}_{codon_number}")) %>%
  select(mutation_name, codon_name)

cor_filt <- cor_df %>% 
  left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
  left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
  filter(codon1 != codon2) %>%
  filter(mut1 %in% parsed_agg$mutation_name,
         mut2 %in% parsed_agg$mutation_name) %>%
  arrange(desc(intra_corr))

obs_df <- fread("results/linkage_out/observed_linkage.gt1000.csv") %>%
  mutate(mut1 = ifelse(mutation1 > mutation2, mutation1, mutation2),
         mut2 = ifelse(mutation1 > mutation2, mutation2, mutation1)) %>%
  select(mut1, mut2, n_obs_comparisons = n_comparisons, obs_corr = corr)

plot_df <- cor_filt %>%
  filter(mut1 %in% parsed_agg$mutation_name,
         mut2 %in% parsed_agg$mutation_name) %>%
  left_join(obs_df) %>%
  filter(!is.na(obs_corr)) %>%
  left_join(parsed_agg %>% select(mut1 = mutation_name, n_mut1 = global_n, max_mut1 = max_prop)) %>%
  left_join(parsed_agg %>% select(mut2 = mutation_name, n_mut2 = global_n, max_mut2 = max_prop)) %>%
  mutate(n_avg = (n_mut1 + n_mut2) / 2)
  # filter(intra_corr > 0)

nrow(plot_df)  
# lr <- lm(obs_corr ~ n_comparisons + n_obs_comparisons + (n_biosamples1 / n_biosamples2) + intra_corr,
#    data = plot_df)

# lr <- lm(obs_corr ~ n_avg + intra_corr,
#          data = plot_df)
# 
# summary(lr)

corr <- cor.test(plot_df$intra_corr, plot_df$obs_corr)
r <- signif(corr$estimate, 2)
p <- signif(corr$p.value, 2)

plot_df %>%
  # filter(n_mut1 > 10000, n_mut2 > 10000) %>%
  ggplot(aes(x = intra_corr, y = obs_corr)) +
  geom_smooth(method = "lm", color = "black", fill = "khaki3") +
  geom_point(color = "black",
             fill = "palegreen4",
             pch = 21,
             size = 3,
             alpha = 0.5) +
  theme_bw() +
  # geom_text_repel(aes(label = str_glue("{mut1}-{mut2}")),
  #           size = 3) +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  labs(x = "Intrahost correlation", y = "Consensus correlation") +
  annotate("text", x = 0, y = 0.9, label = str_glue("r={r},p={p}"),
           hjust = 0)
  # xlim(-1, 1) +
  # ylim(-1, 1)

ggsave("results/linkage_out/intra_versus_obs_scatter.gt1000.pdf", width = 4, height = 4, dpi = 600)

plot_filt <- plot_df %>%
  filter(n_mut1 > 100000, n_mut2 > 100000)

corr <- cor.test(plot_filt$intra_corr, plot_filt$obs_corr)
r <- signif(corr$estimate, 2)
p <- signif(corr$p.value, 2)

plot_filt %>%
  ggplot(aes(x = intra_corr, y = obs_corr)) +
  geom_smooth(method = "lm", color = "black", fill = "khaki3") +
  geom_point(color = "black",
             fill = "palegreen4",
             pch = 21,
             size = 3,
             alpha = 0.5) +
  theme_bw() +
  # geom_text_repel(aes(label = str_glue("{mut1}-{mut2}")),
  #           size = 3) +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  labs(x = "Intrahost correlation", y = "Consensus correlation") +
  annotate("text", x = 0, y = 0.9, label = str_glue("r={r},p={p}"),
           hjust = 0)
# xlim(-1, 1) +
# ylim(-1, 1)

ggsave("results/linkage_out/intra_versus_obs_scatter.gt100000.pdf", width = 4, height = 4, dpi = 600)


plot_df %>% 
  filter(intra_corr > 0.9) %>%
  select(mut1, mut2, intra_corr, obs_corr) %>%
  fwrite("results/linkage_out/candidate_linkage.csv")
# 
lr <- lm(obs_corr ~ log10(n_avg) + intra_corr,
   data = plot_df %>%
     filter(intra_corr > 0) %>%
     filter(n_mut1 > 10000, n_mut2 > 10000))
summary(lr)
relaimpo::calc.relimp(lr)
anova(lr)
# plot_df %>%
#   ggplot(aes(obs_corr, log10(n_avg))) +
#   geom_point() +
#   geom_smooth()
# cor_filt %>%
#   left_join(obs_df) %>%
#   left_join(parsed_agg %>% select(mut1 = mutation_name, n_mut1 = global_n)) %>%
#   left_join(parsed_agg %>% select(mut2 = mutation_name, n_mut2 = global_n)) %>%
#   mutate(n_avg = (n_mut1 + n_mut2) / 2) %>%
#   arrange(desc(intra_corr)) %>% View()

morsels <- foreach(t = seq(1, 100) * 1000) %do% {
  plot_filt <- plot_df %>%
    filter(n_mut1 > t, n_mut2 > t)
  
  plot_filt2 <- plot_filt %>%
    filter(intra_corr > 0)
  
  r_link <- cor(plot_filt$intra_corr, plot_filt$obs_corr)
  r_pos <- cor(plot_filt2$intra_corr, plot_filt2$obs_corr)
  
  tibble(r_link = r_link, r_pos = r_pos, t = t, n_points = nrow(plot_filt))
}  

bind_rows(morsels) %>%
  pivot_longer(!c(n_points, t), names_to = "type", values_to = "r") %>%
  filter(type == "r_link") %>%
  ggplot(aes(x = t, y = r, size = n_points)) +
  geom_point(pch = 21, fill = "indianred",
             color = "black") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold")) +
  labs(x = "No. sequences", y = "r (intrahost vs. consensus)",
       size = "No. mutation pairs")

ggsave("results/linkage_out/r_versus_global_n.pdf", width = 4, height = 4, dpi = 600)

plot_filt <- plot_df %>%
  summarise(n1 = sum(n_mut1 > 1000))
  filter(n_mut1 > 100000, n_mut2 > 100000)

plot_filt %>%
  ggplot(aes(x = intra_corr, y = obs_corr)) +
  geom_smooth(method = "lm", color = "black", fill = "khaki3") +
  geom_point(color = "black",
             fill = "palegreen4",
             pch = 21,
             size = 3,
             alpha = 0.5) +
  theme_bw()
  
cor_df %>% 
  left_join(hookup_df %>% select(mut1 = mutation_name, codon1 = codon_name)) %>%
  left_join(hookup_df %>% select(mut2 = mutation_name, codon2 = codon_name)) %>%
  filter(codon1 != codon2) %>%
  filter(!grepl("stop", mut1)) %>%
  filter(!grepl("stop", mut2)) %>%
  left_join(parsed_agg %>% select(mut1 = mutation_name, n_mut1 = global_n, max_mut1 = max_prop)) %>%
  left_join(parsed_agg %>% select(mut2 = mutation_name, n_mut2 = global_n, max_mut2 = max_prop)) %>%
  mutate(n_mut1 = replace_na(n_mut1, 0)) %>%
  mutate(n_mut2 = replace_na(n_mut2, 0)) %>%
  filter(intra_corr > 0.5|intra_corr < -0.5) %>% 
  summarise(n_total = n(),
            n1 = sum(n_mut1 < 1000 & n_mut2 < 1000),
            n2 = sum(n_mut1 < 1000 & n_mut2 < 1000 & intra_corr > 0),
            n3 = sum(n_mut1 < 1000 & n_mut2 < 1000 & intra_corr > 0 & mut1_max <= 0.5 & mut2_max <= 0.5))
