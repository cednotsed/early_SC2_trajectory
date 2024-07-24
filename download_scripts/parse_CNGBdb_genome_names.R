rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

df <- fread("data/metadata/CNGBdb/CNGBDb_VirusDIP.210324.filt.csv")

# Parse genome names
fna <- readDNAStringSet("data/genomes/CNGBdb/all_CNGBdb.fna")

names(fna) <- gsub("2019-nCoV_", "", names(fna))

to_parse <- names(fna)[!(names(fna) %in% df$virus_name)]

for(acc in to_parse) {
  print(acc)
  
  names(fna)[names(fna) == acc] <- df$virus_name[grepl(acc, df$virus_name)]
}

hookup <- tibble(virus_name = names(fna)) %>%
  left_join(df)

names(fna) <- hookup$accession

writeXStringSet(fna, "data/genomes/all_CNGBdb.renamed.fna")
