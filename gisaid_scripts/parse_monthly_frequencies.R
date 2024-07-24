rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)

file_dir <- "results/allele_frequency_out/monthly_frequencies/"
file_list <- list.files(file_dir, full.names = T)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)
}

bind_rows(morsels) %>%
  fwrite("results/allele_frequency_out/all_monthly_frequencies.csv")
