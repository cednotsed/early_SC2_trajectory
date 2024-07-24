rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

args <- commandArgs(trailingOnly = T)
file_dir <- "data/metadata/gisaid/monthly_mutations/"
file_list <- list.files(file_dir, full.names = T)
file_name <- file_list[args[1]]

id <- gsub(file_dir, "", file_name)
col_month <- str_split(id, "\\.")[[1]][2]
date_chunk <- fread(file_name, sep = "\t") %>%
  mutate(aa_substitutions = gsub("\\(|\\)", "", aa_substitutions))

n_total <- nrow(date_chunk)

# Get master list
aa_list <- paste0(date_chunk$aa_substitutions, collapse = ",")
aa_list <- unique(str_split(aa_list, ",")[[1]])

# Create template matrix
mat <- matrix(0, nrow(date_chunk), length(aa_list))
colnames(mat) <- aa_list

for(i in seq(nrow(date_chunk))) {
  if(i %% 100 == 0){
    print(str_glue("{i}/{n_total}"))
  }
  
  aa_string <- date_chunk[i, ]$aa_substitutions
  
  if(aa_string != "") {
    aa_temp <- str_split(aa_string, ",")[[1]]
    mat[i, aa_temp] <- 1
  }
}

col_sums <- colSums(mat)

temp <- tibble(mutation_name = names(col_sums), 
               n_present = col_sums) %>%
  mutate(n_total = n_total) %>%
  mutate(collection_month = col_month) %>%
  mutate(prop = n_present / n_total)

fwrite(temp, str_glue("results/allele_frequency_out/monthly_frequencies.V2/freq.{col_month}.csv"))

