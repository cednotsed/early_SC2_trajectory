rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

ncbi <- fread("data/metadata/ncbi_virus_150324.complete.excl_lab_env_provirus_vax.gt25000.1Mar20.pango_lineage.csv")
gisaid <- fread("data/metadata/gisaid_150324.1Mar20.pango_lineage.csv") %>%
  separate(taxon, c(NA, "taxon"), "\\|")

merged <- bind_rows(ncbi, gisaid)

fwrite(merged, "data/metadata/all_pango_assignments.csv")
