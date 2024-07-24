rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

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
  separate(lineage, c("main_lineage"), "\\.")

cndb_filt %>%
  ggplot(aes(x = as.Date(sample_collection_date))) +
  geom_histogram(aes(fill = main_lineage),
                 bins = 100) +
  facet_grid(rows = "main_lineage") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", limits = as.Date(c('2019-12-01','2020-04-31'))) +
  labs(x = "Collection date", y = "Genomes", fill = "Lineage",
       title = str_glue("CNGBdb genomes up to 1 Apr 2020 (n={nrow(cndb_filt)})")) 


