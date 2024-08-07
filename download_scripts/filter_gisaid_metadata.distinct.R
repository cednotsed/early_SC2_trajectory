rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

gisaid_meta <- fread("../wuhu_rooting/data/metadata/gisaid/gisaid_metatadata.080724.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x)))

gisaid_filt <- gisaid_meta %>%
  mutate(collection_date = as.Date(collection_date, "%Y-%m-%d")) %>%
  filter(!is.na(collection_date)) %>%
  filter(is.na(`is_low_coverage?`)) %>%
  filter(`is_complete?`) %>%
  filter(!grepl("Manis|Rhinolo", host)) %>%
  mutate(to_keep = grepl("orig|ginal", `passage_details/history`, ignore.case = T) & 
           !grepl("Vero", `passage_details/history`, ignore.case = T)) %>%
  mutate(collection_month = format(collection_date, "%Y-%m")) %>%
  distinct(collection_date, location, .keep_all = T)

fwrite(gisaid_filt, "data/metadata/gisaid/gisaid_metatadata.080724.filt.distinct.tsv")

months <- unique(gisaid_filt$collection_month)

foreach(month = months) %do% {
  temp <- gisaid_filt %>%
    filter(collection_month == month) %>%
    select(aa_substitutions)
  
  print(nrow(temp))
  
  temp %>%
    fwrite(str_glue("data/metadata/gisaid/monthly_mutations.distinct/aa_subs.{month}.080724.tsv"),
           sep = "\t")
}
