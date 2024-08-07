rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

non_nanopore <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv") %>%
  filter(platforms != "OXFORD_NANOPORE")

file_dir <- "results/pipeline_out/parsed_vcfs/"

file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  # file_name <- file_list[10]
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".bcftools.parsed.tsv", "", id)
  
  print(file_name)
  
  temp <- fread(file_name) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(mq = as.numeric(mq)) %>%
    mutate(total_depth = as.numeric(total_depth)) %>%
    select(-qual)
  
  temp %>%
    filter(!is.na(mq)) %>%
    filter(mq >= 30) %>%
    filter(total_depth >= 20) %>%
    filter(alt != "<*>") %>%
    filter(indel_flag == FALSE) %>%
    mutate(id = id, .before = 1)
}

merged <- bind_rows(morsels) %>%
  separate(ad, c("a1", "a2", "a3", "a4"), ",") %>%
  mutate(a1 = replace_na(as.numeric(a1), 0),
         a2 = replace_na(as.numeric(a2), 0),
         a3 = replace_na(as.numeric(a3), 0),
         a4 = replace_na(as.numeric(a4), 0)) %>%
  mutate(p1 = a1 / total_depth,
         p2 = a2 / total_depth,
         p3 = a3 / total_depth,
         p4 = a4 / total_depth)

merged2 <- merged %>% 
  mutate(is_nanopore = id %in% nanopore$biosample) %>%
  select(id, is_nanopore, pos, total_depth, p_A, p_T, p_G, p_C) %>%
  pivot_longer(!c(id, is_nanopore, pos, total_depth), names_to = "allele", values_to = "freq") %>%
  mutate(freq = as.numeric(freq)) %>%
  mutate(allele_depth = total_depth * freq) %>%
  filter(is_nanopore & allele_depth >= 100 | !is_nanopore & allele_depth >= 10)

merged2 %>%
  fwrite("results/allele_frequency_out/merged_frequencies.csv")


  
  
