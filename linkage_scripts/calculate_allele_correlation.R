rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

# fna <- readDNAStringSet("data/alignments/reassembled/reassembled.masked.aln")
# meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")
# non_nanopore <- meta %>%
#   filter(platforms != "OXFORD_NANOPORE")
# 
hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv") %>%
  # filter(!grepl("element|loop", protein_name)) %>%
  select(pos = nucleotide_pos, ref_nuc, var_nuc,
         region, protein_name, ref_AA,
         codon_number, codon_from_gene_start, var_AA,
         mutation_type, mutation_name)
# 
# hookup_select <- hookup %>%
#   filter(mutation_type == "NS") %>%
#   mutate(protein_name = gsub(" \\(part 2\\)| \\(part 1\\)", "", protein_name)) %>%
#   mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
#   mutate(var_AA = ifelse(var_AA == "*", "stop", var_AA)) %>%
#   mutate(protein_name = case_when(protein_name == "S" ~ "Spike",
#                                   protein_name %in% c("NS3a", "NSP3a", "NSP12a", "NS8a") ~ gsub("a", "", protein_name),
#                                   T ~ protein_name)) %>%
#   mutate(mutation_name = str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}")) %>%
#   filter(!grepl("\\*", ref_AA)) %>% # Remove stop codon mutations 
#   distinct(mutation_name, protein_name, codon_number)

df <- fread("results/allele_frequency_out/merged_frequencies.csv") %>%
  filter(freq > 0.1)

parsed <- df %>% 
  select(id, pos, total_depth, allele, freq) %>%
  mutate(var_nuc = gsub("p_", "", allele)) %>%
  left_join(hookup) %>%
  filter(mutation_type == "NS") %>%
  mutate(mutation_name = str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}")) %>%
  filter(!grepl("\\*", ref_AA)) # Remove stop codon mutations

mat <- parsed %>%
  group_by(id, mutation_name) %>%
  summarise(sum_freq = sum(freq)) %>%
  ungroup() %>%
  pivot_wider(id_cols = id, names_from = mutation_name, values_from = sum_freq) %>%
  column_to_rownames("id")

mat[is.na(mat)] <- 0

# Get mutation counts
col_sums <- colSums(mat > 0)
count_df <- tibble(mutation = names(col_sums), 
                   n_biosamples = col_sums)
to_keep <- names(col_sums)[col_sums > 10]

mat_filt <- mat[, to_keep]

# cor_mat <- cor(mat_filt)
# 
# cor_mat[upper.tri(cor_mat)] <- NA
# 
# rownames(cor_mat) <- colnames(mat_filt)
# 
# cor_df <- as.data.frame(cor_mat) %>%
#   rownames_to_column("mutation1") %>%
#   pivot_longer(!mutation1, names_to = "mutation2", values_to = "corr") %>%
#   filter(mutation1 != mutation2) %>%
#   filter(!is.na(corr)) %>%
#   arrange(desc(corr))

pairs <- combn(to_keep, 2)

print(ncol(pairs))

morsels <- foreach(i = seq(ncol(pairs))) %do% {
  print(i)
  pair <- pairs[, i]
  pair <- c("Spike_D614G", "NSP12_P323L")
  mat_temp <- mat[, pair]
  samples_to_keep <- rowSums(mat_temp) != 0
  mat_temp <- mat_temp[samples_to_keep, ]
  corr <- cor(mat_temp[, 1], mat_temp[, 2])
  tibble(mutation1 = pair[1], mutation2 = pair[2], corr = corr)
}

merged <- bind_rows(morsels) %>%
  arrange(desc(corr)) %>%
  left_join(hookup %>% dplyr::rename(mutation1 = mutation_name,
                                            protein_name1 = protein_name,
                                            codon_number1 = codon_number)) %>%
  left_join(hookup %>% dplyr::rename(mutation2 = mutation_name,
                                            protein_name2 = protein_name,
                                            codon_number2 = codon_number)) %>%
  
  filter(!(protein_name1 == protein_name2 & codon_number1 == codon_number2)) %>%
  mutate(rsquared = corr * corr) %>%
  left_join(count_df %>% dplyr::rename(mutation1 = mutation, n_biosamples1 = n_biosamples)) %>%
  left_join(count_df %>% dplyr::rename(mutation2 = mutation, n_biosamples2 = n_biosamples))

merged_filt <- merged %>%
  filter(rsquared > 0.5)

fwrite(merged, "results/allele_frequency_out/linkage.csv")
fwrite(merged_filt, "results/allele_frequency_out/linkage.r50.csv")

# gisaid_meta <- fread("data/metadata/gisaid/gisaid_metatadata.080724.filt.tsv")
# 
# gisaid_filt <- gisaid_meta %>%
#   filter(!(collection_month %in% c("2019-12", "2020-01", "2020-02", "2020-03")))
# 
# res_morsels <- foreach(i = seq(nrow(merged_filt))) %do% {
#   print(i)
#   row <- merged_filt[i, ]
#   mut1 <- row$mutation1
#   mut2 <- row$mutation2
#   
#   res <- gisaid_filt %>%
#     summarise(sum_linked = sum(grepl(mut1, aa_substitutions) & grepl(mut2, aa_substitutions)),
#               sum_unlinked = sum((grepl(mut1, aa_substitutions) & !grepl(mut2, aa_substitutions)) | 
#                                    (grepl(mut2, aa_substitutions) & !grepl(mut1, aa_substitutions))))
# 
#   temp <- row %>%
#     mutate(sum_linked = res$sum_linked,
#            sum_unlinked = res$sum_unlinked)
#   temp %>%
#     mutate(ratio = sum_linked / sum_unlinked) %>%
#     select(ratio)
#   return(temp)
# }

