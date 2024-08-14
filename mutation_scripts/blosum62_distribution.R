rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
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

# Add blosum score
data(BLOSUM62)

changes <- parsed_agg %>%
  distinct(ref_AA, var_AA)

# changes %>%
  # filter(ref_AA == "i")
change_morsels <- foreach(i = seq(nrow(changes))) %do% {
  # i = 218
  row <- changes[i, ]
  
  row %>%
    mutate(blosum62_score = BLOSUM62[row$ref_AA, ifelse(row$var_AA == "stop", "*", row$var_AA)])
}

change_df <- bind_rows(change_morsels)

plot_df <- parsed_agg %>%
  mutate(is_fixed = max_prop > 0.9) %>%
  filter(ref_AA %in% possible_variants,
         var_AA %in% possible_variants) %>%
  left_join(change_df)

plot_df %>%
  ggplot(aes(x = factor(blosum62_score), y = log10(global_n), fill = is_fixed)) +
  geom_boxplot() +
  facet_grid(rows = vars(is_fixed),
             scales = "free") +
  labs(x = "BLOSUM62 score", y = "Log10(no. genomes)") +
  scale_fill_manual(values = c("indianred", "turquoise3")) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))

ggsave("results/mutation_out/global_n_versus_blosum62.pdf", dpi = 600, width = 4, height = 4)

linreg <- lm(log10(global_n) ~ blosum62_score,
   data = plot_df %>% filter(!is_fixed))

test_df <- plot_df %>%
  filter(is_fixed)

cor.test(test_df$global_n, test_df$blosum62_score, method = "spearman")

hist(linreg$residuals)
summary(linreg)

plot_df %>%
  ggplot(aes(x = blosum62_score, y = log10(global_n), color = is_fixed)) +
  geom_smooth(method = "lm") +
  geom_point()

# mat <- monthly_agg %>% 
#   arrange(collection_month) %>%
#   pivot_wider(id_cols = mutation_name, names_from = collection_month, values_from = monthly_prop) %>%
#   column_to_rownames("mutation_name")
  
# mat[is.na(mat)] <- 0
# 
# cor(t(mat))
plot_df %>%
  group_by(blosum62_score) %>%
  summarise(n = n_distinct(mutation_name)) %>%
  ggplot(aes(x = factor(blosum62_score), y = n, fill = factor(blosum62_score))) +
  geom_bar(stat = "identity",
           color = "black") +
  scale_fill_viridis_d() +
  labs(x = "BLOSUM62 score", y = "Distinct mutations observed") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"))
  
ggsave("results/mutation_out/blosum_counts.pdf", dpi = 600, width = 4, height = 4)

# as.data.frame(blosum62) %>%
#   rownames_to_column("var1") %>%
#   pivot_longer(!var1, names_to = "var2", values_to = "score") %>%
#   filter(var1 != "")
#   arrange(desc(score)) %>% 

hookup %>%
  distinct(region, protein_name, codon_number) %>%
  filter(codon_number != -1) %>% 
  group_by(region) %>%
  summarise(n = n()) %>% View()
plot_df %>%
  mutate(ORF = ifelse(grepl("NSP", protein_name), "ORF1ab", protein_name)) %>%
  group_by(ORF) %>%
  summarise(n = n()) %>%
  View()
  
hookup %>% 
  filter(protein_name == "Spike") %>% View()

plot_df %>%
  mutate(conservative = blosum62_score >= 0) %>%
  group_by(is_fixed, conservative) %>%
  summarise(n = n()) %>%
  group_by(is_fixed) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = is_fixed, y = prop, fill = conservative)) +
  geom_bar(stat = "identity")


mat <- plot_df %>%
  mutate(conservative = blosum62_score >= 0) %>%
  group_by(is_fixed, conservative) %>%
  summarise(n = n()) %>%
  pivot_wider(id_cols = is_fixed, names_from = conservative, values_from = n) %>%
  column_to_rownames("is_fixed")

fisher.test(mat)

