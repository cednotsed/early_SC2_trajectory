rm(list = ls())
setwd("/SAN/ugi/HAP_VAP/early_SC2_trajectory")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

meta <- fread("data/metadata/all_sra_metadata.csv")
aa_meta <- fread("data/metadata/gisaid/all_aa_freq.080724.csv") %>%
  select(collection_month)

parsed_agg <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")
parsed_filt <- parsed_agg %>%
  filter(global_n > 1000)

savs <- parsed_filt$mutation_name

file_dir <- "data/metadata/gisaid/presence_matrix.temp/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list[1:500]) %do% {
  # file_name = file_list[1]
  mut <- gsub(file_dir, "", file_name)
  mut <- gsub(".presence.csv.gz", "", mut)
  
  print(mut)
  temp <- fread(file_name)
  colnames(temp) <- mut
  gc()
  return(temp)
}

mat <- bind_cols(aa_meta, morsels)

fwrite(mat, "data/metadata/gisaid/presence_absence_matrix.csv")
mat[1:5, 1:5]


