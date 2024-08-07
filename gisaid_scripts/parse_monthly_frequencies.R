rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

file_dir <- "results/allele_frequency_out/monthly_frequencies.V2.distinct/"
# file_dir <- "results/allele_frequency_out/monthly_frequencies.V2/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name) %>%
    mutate(across(everything(), as.character))
}

bind_rows(morsels) %>% 
  filter(mutation_name != "") %>%
  fwrite("results/allele_frequency_out/all_monthly_frequencies.distinct.csv")
  # fwrite("results/allele_frequency_out/all_monthly_frequencies.csv")
