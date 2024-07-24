rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.csv") %>%
  filter(!grepl("element|loop", protein_name)) %>%
  mutate(protein_name = gsub(" \\(part 2\\)| \\(part 1\\)", "", protein_name)) %>%
  mutate(protein_name = gsub("ORF", "NS", protein_name)) %>%
  mutate(var_AA = ifelse(var_AA == "*", "stop", var_AA)) %>%
  mutate(protein_name = case_when(protein_name == "S" ~ "Spike",
                                  protein_name %in% c("NS3a", "NSP3a", "NSP5a", 
                                                      "NSP12a", "NSP15a", "NS8a") ~ gsub("a", "", protein_name),
                                  T ~ protein_name)) %>%
  mutate(mutation_name = ifelse(mutation_type == "NS", 
                                str_glue("{protein_name}_{ref_AA}{codon_number}{var_AA}"),
                                str_glue("{protein_name}_{ref_nuc}{nucleotide_pos}{var_nuc}")))
         
fwrite(hookup, "data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")
