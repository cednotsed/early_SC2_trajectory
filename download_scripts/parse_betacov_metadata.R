rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

gisaid_meta <- fread("data/metadata/betacoronavirus/gisaid_betacoronavirus.complete.dates_and_loc.040424.tsv") %>%
  left_join(fread("data/metadata/betacoronavirus/gisaid_betacoronavirus.complete.sequencing_tech.040424.tsv")) %>%
  left_join(fread("data/metadata/betacoronavirus/gisaid_betacoronavirus.complete.patient.040424.tsv")) %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(database = "GISAID") %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
  separate(location, c(NA, "country"), " / ", remove = F) %>%
  select(database, accession = accession_id, collection_date,
         submission_date, virus_name, location, 
         host, isolation_source = specimen, sequencing_technology) %>%
  mutate(genus = "Betacoronavirus")

ncbi <- fread("data/metadata/betacoronavirus/coronaviridae.taxid1118.gt25000.excl_cov2_provirus_env_lab_host_vax.020424.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(database = "NCBI") %>%
  select(database, accession, collection_date,
         virus_name = genbank_title, location = geo_location, 
         country, submitting_lab = organization, host, isolation_source,
         nuc_completeness, sra_accession, biosample, bioproject, genus)

merged <- gisaid_meta %>%
  bind_rows(ncbi)

fwrite(merged, "data/metadata/betacoronavirus/all_coronaviruses.csv")

  

