rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)

merged <- fread("results/allele_frequency_out/merged_frequencies.csv")
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>%
  filter(!grepl("element|loop", protein_name))

fna <- readDNAStringSet("data/alignments/reassembled/reassembled.masked.aln")
meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")

non_nanopore <- meta %>%
  filter(platforms != "OXFORD_NANOPORE")

nanopore <- meta %>%
  filter(platforms == "OXFORD_NANOPORE")

merged %>%
  filter(id %in% names(fna)) %>%
  filter(id %in% non_nanopore$biosample) %>%
  summarise(n = n_distinct(id))

allele_count <- merged %>%
  filter(id %in% non_nanopore$biosample) %>%
  filter(id %in% names(fna)) %>%
  filter(freq > 0.1 & freq < 0.5) %>%
  group_by(pos, allele) %>%
  summarise(n = n_distinct(id)) %>%
  mutate(var_nuc = gsub("p_", "", allele)) %>%
  arrange(desc(n))

multi_alleles <- allele_count %>%
  filter(n > 1) %>%
  mutate(pos = as.numeric(pos)) %>%
  left_join(hookup %>%
              select(pos = nucleotide_pos, ref_nuc, var_nuc, 
                     region, protein_name, ref_AA, 
                     codon_number, codon_from_gene_start, var_AA, 
                     mutation_type, nucleotide_length)) %>%
  mutate(mutation_name = ifelse(mutation_type == "NS", 
                                str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"),
                                str_glue("{protein_name}_{ref_nuc}{pos}{var_nuc}"))) %>%
  mutate(is_reference = ifelse(ref_nuc == var_nuc, "Reference","Non-reference")) %>%
  mutate(mutation_name = ifelse(n > 10, mutation_name, NA)) %>%
  filter(n > 1) %>%
  mutate(protein_name = factor(protein_name, unique(hookup$protein_name)))

pal <- distinctColorPalette(n_distinct(multi_alleles$region))
pal <- c("#DCAC51", "#CB9E88", "#B44CE7", "#DEE561", "#DDAAD2", "#D5DCD9", "#D1E0A2", "#897FD4", "#84B2D8",
         "#7ADD91", "#DC6AC2", "#90EA4F", "#E45B63", "#70D7D0")
pal <- setNames(pal, unique(multi_alleles$region))

multi_alleles %>%
  ggplot(aes(x = pos, y = n, fill = region, shape = mutation_type)) +
  geom_point(color = "black", 
             alpha = 0.8) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = pal, guide = "none") +
  geom_text_repel(aes(label = mutation_name)) +
  facet_grid(rows = vars(is_reference)) +
  labs(x = "Position on Wuhan-Hu-1", y = "No. biosamples", shape = "Mutation type") +
  theme_bw() 

multi_alleles %>%
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
