rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(ggrepel)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  distinct(mutation_name, protein_name, codon_number, region, ref_AA, var_AA)

merged <- fread("results/allele_frequency_out/all_sites.significant.csv")
parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv") %>%
  dplyr::select(mutation_name, max_prop)

plot_df <- merged %>%
  left_join(parsed_agg) %>%
  left_join(hookup) %>%
  filter(!is.na(max_prop)) %>%
  mutate(fixed_annot = ifelse(max_prop > 0.1, str_glue("{ref_AA}{codon_number}{var_AA}"), NA),
         is_fixed = max_prop > 0.9) %>%
  mutate(fixed_annot = ifelse(region == "ORF1ab", str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"), fixed_annot))

# pal <- c("#DCAC51", "#CB9E88", "#B44CE7", "#DEE561", "#DDAAD2", "#D5DCD9", "#D1E0A2", "#897FD4", "#84B2D8",
         "#7ADD91", "#DC6AC2", "#90EA4F", "#E45B63", "#70D7D0")
# pal <- setNames(pal, unique(plot_df$region))
pal <- list(ORF1ab = "#84B2D8", S = "#7ADD91", ORF3a = "#DEE561", 
            M = "#CB9E88", ORF7a = "#D5DCD9", ORF8 = "#897FD4", 
            N = "#B44CE7")

plot_df %>%
  mutate(region = factor(region, c("ORF1ab", "S", "ORF3a", "M", "ORF7a", "ORF8", "N"))) %>%
  filter(max_prop > 0.1) %>%
  ggplot(aes(x = region, y = max_prop)) +
  geom_point(aes(fill = region, shape = is_fixed),
             position = position_jitter(seed = 66),
             color = "black", 
             alpha = 0.8,
             size = 3) +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_manual(values = pal) +
  geom_text_repel(aes(label = fixed_annot), 
                  position = position_jitter(seed = 66),
                  size = 4,
                  ylim = c(0.2, 0.9)) +
  theme_bw() +
  theme(legend.position = "none", 
        # text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Position on Wuhan-Hu-1", y = "No. biosamples detected", shape = "Max. monthly freq.")

ggsave("results/allele_frequency_out/mutations_detected.pdf", dpi = 600, height = 3, width = 12)

plot_df %>%
  group_by(protein_name, region, nucleotide_length) %>%
  summarise(n = n_distinct(mutation_name)) %>%
  filter(nucleotide_length > 100) %>%
  mutate(protein_name = factor(protein_name, unique(hookup$protein_name))) %>%
  mutate(normalised_count = n / nucleotide_length) %>%
  ggplot(aes(x = protein_name, y = normalised_count, fill = region)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Region", y = "Mutations / region length", fill = "ORF")

ggsave("results/allele_frequency_out/mutations_detected.pdf", dpi = 600, height = 5, width = 12)
# hookup %>%
#   group_by(protein_name) %>%
#   summarise(min= min(nucleotide_pos),
#             max = max(nucleotide_pos)) %>% View()

dnds_meta <- hookup %>%
  filter(ref_nuc != var_nuc) %>%
  filter(pos_in_codon != -1) %>%
  group_by(nucleotide_pos, protein_name) %>%
  summarise(n_ns = any(mutation_type == "NS")) %>%
  group_by(protein_name) %>%
  summarise(total_ns = sum(n_ns),
            total_s = sum(!n_ns))

multi_alleles %>%
  group_by(protein_name, region) %>%
  summarise(n_ns = n_distinct(pos[mutation_type == "NS"]),
            n_s = n_distinct(pos[mutation_type == "S"])) %>%
  left_join(dnds_meta) %>%
  mutate(pnps = (n_ns / total_ns) / (n_s / total_s)) %>%
  filter(!is.na(total_ns),
         pnps != Inf) %>% 
  mutate(protein_name = factor(protein_name, unique(hookup$protein_name))) %>%
  ggplot(aes(x = protein_name, y = pnps, fill = region)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Region", y = "pn/ps", fill = "ORF")
