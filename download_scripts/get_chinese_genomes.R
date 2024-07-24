rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

cndb <- fread("data/metadata/CNGBdb/metadata.220324.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(data_source %in% c("NMDC", "Genome Warehouse", "CNGBdb", "GenBase"))

cndb %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  distinct(data_source)

cndb_filt <- cndb %>%
  filter(related_id == "null") %>% 
  filter(sequence_length > 29000) %>% 
  filter(sample_collection_date <=  as.Date("2020-03-01")) %>% 
  select(database = data_source, accession = accession_id, virus_name = virus_strain_name,
         location = location,
         cncb_completeness = nuc.completeness, cncb_quality = sequence_quality, cncb_quality_string = quality_assessment)

fna <- c(readDNAStringSet("data/genomes/chinese_genomes/CNGB.fasta"),
         readDNAStringSet("data/genomes/chinese_genomes/GWH.fasta"),
         readDNAStringSet("data/genomes/chinese_genomes/NMDC.fasta"))

cndb_filt %>%
  filter(!(accession %in% names(fna)))
