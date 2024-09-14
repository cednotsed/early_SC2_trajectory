rm(list = ls())
setwd("/SAN/ugi/HAP_VAP/early_SC2_trajectory")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)

parsed_agg <- fread("results/allele_frequency_out/observed_out/all_time.aggregate.csv")
parsed_filt <- parsed_agg %>%
  filter(global_n > 1000)

savs <- parsed_filt$mutation_name

aa_meta <- fread("data/metadata/gisaid/all_aa_freq.080724.csv")
aa_strings <- aa_meta$aa_substitutions

print(length(savs))

cl <- makeCluster(24)
registerDoParallel(cl)

foreach(mut = savs, .packages = c("tidyverse", "data.table")) %dopar% {
  # mut = savs[1]
  print(mut)
  presence <- grepl(mut, aa_strings)
  
  tibble(mut = presence) %>%
    fwrite(str_glue("data/metadata/gisaid/presence_matrix.temp/{mut}.presence.csv.gz"))
  gc()

  return(NULL)
}

stopCluster(cl)

