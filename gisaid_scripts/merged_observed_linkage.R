rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

file_dir <- "results/linkage_out/observed_linkage.early.temp/"

file_list <- list.files(file_dir, full.names = T)

length(file_list)

morsels <- foreach(file_name = file_list) %do% {
  fread(file_name)  
}

bind_rows(morsels) %>%
  fwrite("results/linkage_out/observed_linkage.early.gt100.csv")
