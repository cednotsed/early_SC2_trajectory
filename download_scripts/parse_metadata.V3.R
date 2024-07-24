rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

# Chinese genomes
cncb <- fread("data/metadata/CNCB/metadata.220324.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(data_source %in% c("NMDC", "Genome Warehouse", "CNGBdb", "GenBase")) %>%
  filter(related_id == "null")

gisaid <- fread("data/metadata/gisaid/gisaid.dates_and_loc.complete.1Mar20.040424.tsv") %>%
  inner_join(fread("data/metadata/gisaid/gisaid.patient.complete.1Mar20.040424.tsv")) %>%
  inner_join(fread("data/metadata/gisaid/gisaid.sequencing_tech.complete.1Mar20.040424.tsv")) %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  rename_all(~gsub("\\?", "", .x)) %>%
  mutate(database = "GISAID") %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
  separate(location, c(NA, "country"), " / ", remove = F)

ncbi <- fread("data/metadata/NCBI/NCBI.1Mar20.gt29000.040424.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
  mutate(database = "NCBI")

# Filter
cncb_filt <- cncb %>%
  filter(sequence_length > 29000) %>% 
  filter(sample_collection_date <=  as.Date("2020-03-01")) %>%
  separate(location, c("country"), " / ", remove = F) %>%
  select(database = data_source, accession = accession_id, 
         collection_date = sample_collection_date, virus_name = virus_strain_name,
         location = location, country, submitting_lab, 
         nuc_completeness = nuc.completeness, 
         cncb_quality = sequence_quality, 
         cncb_quality_string = quality_assessment, 
         sequence_length)

gisaid_filt <- gisaid %>%
  filter(!grepl("Manis|Rhinolo", host)) %>%
  separate(location, c(NA, "country"), " / ", remove = F) %>%
  select(database, accession = accession_id, collection_date,
         submission_date, virus_name, location, country,
         host, isolation_source = specimen, sequencing_technology)

ncbi_filt <- ncbi %>%
  filter(accession != "NC_045512.2") %>%
  select(database, accession, collection_date,
         virus_name = isolate, location = geo_location, 
         country, submitting_lab = organization, host, isolation_source,
         nuc_completeness, sra_accession, biosample, bioproject)

merged <- cncb_filt %>%
  bind_rows(gisaid_filt, ncbi_filt) 

# Remove missing collection month
merged_filt <- merged %>%
  separate(collection_date, c("Y", "M", "D"), "-", remove = F) %>%
  filter(!is.na(M))

fwrite(merged_filt, "data/metadata/all_metadata.040424.tsv", 
       sep = "\t")

# Extract genomes
gisaid_fna <- readDNAStringSet("data/genomes/gisaid/gisaid.complete.1Mar20.040424.fna")
names(gisaid_fna) <- str_split(names(gisaid_fna), "\\|", simplify = T)[, 2]

chinese_fna <- c(readDNAStringSet("data/genomes/CNCB/CNGB.fasta"),
                 readDNAStringSet("data/genomes/CNCB/GWH.fasta"))
names(chinese_fna) <- str_split(names(chinese_fna), "\t", simplify = T)[, 1]

merged_fna <- c(gisaid_fna, chinese_fna,
                readDNAStringSet("data/genomes/NCBI/NCBI.1Mar20.gt29000.040424.fna"),
                readDNAStringSet("data/genomes/CNCB/NMDC.fasta"))

merged %>%
  filter(!(accession %in% names(merged_fna)))

merged_fna_filt <- merged_fna[merged_filt$accession]

length(merged_fna_filt) == nrow(merged_filt)

writeXStringSet(merged_fna_filt, "data/genomes/all_genomes.040424.fna")
