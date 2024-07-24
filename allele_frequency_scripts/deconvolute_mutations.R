rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

fna <- readDNAStringSet("data/alignments/reassembled/reassembled.masked.aln")
meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  select(pos = nucleotide_pos, ref_nuc, var_nuc, 
         region, protein_name, ref_AA, 
         codon_number, pos_in_codon, codon_from_gene_start, var_AA, 
         mutation_type, mutation_name)

non_nanopore <- meta %>%
  filter(platforms != "OXFORD_NANOPORE")

df <- fread("results/allele_frequency_out/merged_frequencies.csv")

# # Check all allele frequency entries are unique
# df %>%
#   filter(id %in% names(fna)) %>%
#   filter(id %in% non_nanopore$biosample) %>%
#   group_by(id, pos, allele) %>%
#   summarise(n = n()) %>%
#   arrange(desc(n))

df_filt <- df %>%
  filter(id %in% names(fna)) %>%
  filter(id %in% non_nanopore$biosample) %>%
  mutate(pos = as.numeric(pos)) %>%
  mutate(var_nuc = gsub("p_", "", allele)) %>%
  left_join(hookup %>% select(ref_nuc, pos, var_nuc, protein_name, pos_in_codon, ref_AA, codon_number))

# Check if there are any double/triple mutations
codon_allele_count <- df_filt %>%
  filter(ref_nuc != var_nuc) %>%
  group_by(protein_name, codon_number) %>%
  group_by(id, protein_name, codon_number, ref_AA) %>%
  summarise(n_codon_pos = n_distinct(pos_in_codon),
            n_alleles = n()) %>%
  arrange(desc(n_codon_pos)) %>% 
  ungroup() 

# Deconvolute multi-position mutations
# Chuck codons that have more than one allele per position
triple_df <- codon_allele_count %>%
  filter(n_codon_pos == 3) %>%
  filter(n_alleles == n_codon_pos)

double_df <- codon_allele_count %>%
  filter(n_codon_pos == 2) %>%
  filter(n_alleles == n_codon_pos)

single_df <- codon_allele_count %>%
  filter(n_codon_pos == 1)

# Parse single-position mutations
parsed_single <- single_df %>%
  left_join(df_filt) %>%
  left_join(hookup %>% select(ref_nuc, pos, var_nuc, 
                              protein_name, ref_AA, codon_number, 
                              var_AA, pos_in_codon, mutation_type)) %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}")) %>%
  filter(ref_nuc != var_nuc,
         mutation_type == "NS") %>%
  group_by(id, mutation_name) %>%
  summarise(freq = sum(freq))

fwrite(parsed_single, "results/allele_frequency_out/single_position_sites.raw.csv")

parsed_single %>%
  group_by(mutation_name) %>%
  summarise(n = n_distinct(id),
            n_filt = n_distinct(id[freq > 0.1]),
            median_freq = median(freq),
            max_freq = max(freq)) %>%
  arrange(desc(n)) %>%
  fwrite("results/allele_frequency_out/allele_frequencies/single_position_sites.csv")

# Treat all mutations as single
parsed <- df_filt %>%
  left_join(hookup) %>%
  filter(mutation_type == "NS",
         ref_nuc != var_nuc) %>%
  group_by(mutation_name) %>%
  group_by(id, mutation_name) %>%
  summarise(freq = sum(freq)) %>%
  group_by(mutation_name) %>%
  summarise(n = n_distinct(id),
            n_filt = n_distinct(id[freq > 0.1]),
            median_freq = median(freq),
            max_freq = max(freq)) %>%
  arrange(desc(n))

fwrite(parsed, "results/allele_frequency_out/allele_frequencies/all_sites.csv")

# # Parse double-position mutations
# codon_df <- hookup %>%
#   distinct(protein_name, ref_AA, ref_nuc, codon_number, pos_in_codon)
# 
# double_morsels <- foreach(i = seq(nrow(double_df))) %do% {
#   print(i)
#   row <- double_df[i, ]
#   
#   temp <- row %>%
#     left_join(df_filt) %>%
#     filter(ref_nuc != var_nuc)
#   
#   codon_temp <- codon_df %>%
#     filter(protein_name == row$protein_name, 
#            codon_number == row$codon_number) %>%
#     arrange(pos_in_codon) %>%
#     left_join(temp %>%
#                 select(pos_in_codon, var_nuc)) %>%
#     mutate(var_nuc = ifelse(is.na(var_nuc), ref_nuc, var_nuc))
#   
#   var_AA <- paste0(codon_temp$var_nuc, collapse = "")
#   var_AA <- as.character(translate(DNAString(var_AA)))
#   freq <- mean(temp$freq)
#   
#   final <- temp %>%
#     distinct(id, protein_name, ref_AA, codon_number) %>%
#     mutate(var_AA = var_AA, freq = freq)
#   
#   return(final)
# }
# 
# # Parse triple-position mutations
# triple_morsels <- foreach(i = seq(nrow(triple_df))) %do% {
#   print(i)
#   row <- triple_df[i, ]
#   temp <- row %>%
#     left_join(df_filt) %>%
#     filter(ref_nuc != var_nuc) %>%
#     arrange(pos_in_codon)
#   
#   # Check that alleles are in normal order
#   if(all(temp$pos_in_codon == c(1, 2, 3))) {
#     var_AA <- paste0(temp$var_nuc, collapse = "")
#     var_AA <- as.character(translate(DNAString(var_AA)))
#     freq <- mean(temp$freq)
#     
#     final <- temp %>%
#       distinct(id, protein_name, ref_AA, codon_number) %>%
#       mutate(var_AA = var_AA, freq = freq)
#     
#     return(final)
#   }
# }
# 
# 
# double_parsed <- bind_rows(double_morsels) %>%
#   mutate()
# triple_parsed <- bind_rows(triple_morsels)
# 
# double_parsed
