rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

ncbi <- fread("data/metadata/archive/ncbi_virus_150324.complete.excl_lab_env_provirus_vax.gt25000.1Mar20.csv")
ncbi %>% View()
meta <- fread("data/metadata/parsed_metadata.ABC.csv")
pekar <- fread("data/metadata/literature/Pekar_science.abp8337_data_s1.tsv") %>%
  select(accession = `Accession ID`, analysis = Analysis) %>%
  bind_rows(fread("data/metadata/literature/Pekar_science.abp8337_data_s2.tsv") %>%
              dplyr::rename(accession = Accession, location = Location, analysis = Analysis)) %>%
  filter(grepl("Early|Diamond|NYC", analysis))

# Before 2 Mar 2020
fj <- readDNAStringSet("data/metadata/literature/FountainJonesVE_suppl_alignment.fasta")

fj_df <- tibble(genome_name = names(fj)) %>% 
  separate(genome_name, c(NA, "accession", "collection_date"), "\\|") %>%
  filter(collection_date < as.Date("2020-03-02"))

shanghai <- fread("data/metadata/literature/shanghai_paper.csv")

existing <- unique(c(shanghai$accession, fj_df$accession, pekar$accession))

existing[!(existing %in% gsub("\\.1", "", meta$accession))]

meta$accession


cndb <- fread("data/metadata/CNGBdb/CNGBdb_VirusDIP.210324.csv")
cndb
View(cndb)

