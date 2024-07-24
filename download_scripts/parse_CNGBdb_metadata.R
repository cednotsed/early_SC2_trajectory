rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

df <- fread("data/metadata/CNGBdb/CNGBdb_VirusDIP.210324.csv")

cov_df <- fread("data/metadata/CNGBdb/CNGBDb_VirusDIP.210324.curated.txt")

df_filt <- df %>% 
  rename_all(~tolower(gsub(" ", "_", .x))) %>% 
  filter(data_source_platform == "CNGB") %>%
  separate(location, c("country", "province"), ": ") %>%
  mutate(province = ifelse(province == "Weifang", "Shandong", province)) %>%
  mutate(province = ifelse(province == "Wuhan", "Hubei", province)) %>%
  left_join(cov_df) %>%
  dplyr::rename(accession = sequence_id,
                collection_date = sample_collection_date) %>%
  filter(collection_date <= "2020-03-01")

fwrite(df_filt, "data/metadata/CNGBdb/CNGBDb_VirusDIP.210324.filt.csv")


