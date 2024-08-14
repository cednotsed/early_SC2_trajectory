rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(randomcoloR)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
annot_df <- fread("data/metadata/wuhan-hu-1_genome_annotations_V2.csv")
protein_lengths <- annot_df %>%
  filter(region_type == "coding") %>% 
  group_by(region) %>%
  summarise(protein_length = sum(protein_length)) %>%
  mutate(protein_length = ifelse(region == "ORF1ab", 7097, protein_length)) %>%
  # Remove stop codons
  mutate(protein_length = protein_length - 1) %>%
  dplyr::rename(region_name = region)

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")

fixed <- parsed_agg %>%
  filter(max_prop > 0.9)

drift <- parsed_agg %>%
  filter(max_prop < 0.1)

major <- parsed_agg %>%
  filter(max_prop > 0.1 & max_prop < 0.9)

nrow(fixed) / nrow(parsed_agg) * 100
nrow(drift) / nrow(parsed_agg) * 100
nrow(major) / nrow(parsed_agg) * 100
# Remove non-canonical amino acids
possible_variants <- c("G", "A", "V", "L", "I",
                       "T", "S", "M", "C", "P",
                       "F", "Y", "W", "H", "K",
                       "R", "D", "E", "N", "Q",
                       "*")

plot_df <- parsed_agg %>%
  mutate(mut_cat = case_when(mutation_name %in% fixed$mutation_name ~ ">90%",
                             mutation_name %in% drift$mutation_name ~ "<10%",
                             mutation_name %in% major$mutation_name ~ "10-90%")) %>%
  filter(ref_AA %in% possible_variants,
         var_AA %in% possible_variants) %>%
  mutate(region_name = ifelse(grepl("NSP", protein_name), "ORF1ab", protein_name)) %>%
  group_by(region_name, mut_cat) %>%
  summarise(n = n()) %>%
  mutate(region_name = gsub("NS", "ORF", region_name)) %>%
  mutate(region_name = ifelse(region_name == "ORF3", "ORF3a", region_name)) %>%
  mutate(region_name = ifelse(region_name == "Spike", "S", region_name)) %>%
  left_join(protein_lengths) %>%
  mutate(region_name = factor(region_name, c("ORF1ab", "S", "ORF3a",
                                             "E", "M", "ORF6", 
                                             "ORF7a", "ORF7b", "ORF8",
                                             "N", "ORF10")))

pal <- c("#DCAC51", "#CB9E88", "#B44CE7", "#DEE561", "#DDAAD2", "#D5DCD9", "#D1E0A2", "#897FD4", "#84B2D8",
         "#7ADD91", "#DC6AC2", "#90EA4F", "#E45B63", "#70D7D0")
pal <- setNames(pal, unique(plot_df$region_name))

plot_df %>%
  mutate(mut_cat = factor(mut_cat, c("<10%", "10-90%", ">90%"))) %>%
  ggplot(aes(x = region_name, y = n / protein_length, fill = region_name)) +
  geom_bar(stat = "identity",
           color = "black") +
  geom_text(aes(label = n),
            angle = 90,
            position = position_stack(vjust = 0.5)) +
  facet_grid(rows = vars(mut_cat),
             scales = "free") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Region", y = "Distinct mutations / protein length")
 
ggsave("results/mutation_out/mutation_by_ORF.pdf", dpi = 600, width = 4, height = 4)

