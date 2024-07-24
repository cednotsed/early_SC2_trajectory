rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)

args <- commandArgs(trailingOnly=TRUE)

df <- fread("results/allele_frequency_out/linkage.r50.csv") %>%
  arrange(desc(rsquared))

gisaid_meta <- fread("data/metadata/gisaid/gisaid_metatadata.080724.filt.tsv")

gisaid_filt <- gisaid_meta %>%
  filter(!(collection_month %in% c("2019-12", "2020-01", "2020-02", "2020-03")))

# gisaid_filt <- gisaid_filt %>%
#   sample_n(100000)

i <- args[1]

row <- df[i, ]
corr <- row$corr

mut1 <- row$mutation1
mut2 <- row$mutation2

res <- gisaid_filt %>%
  summarise(sum_linked = sum(grepl(mut1, aa_substitutions) & grepl(mut2, aa_substitutions)),
            sum_unlinked = sum((grepl(mut1, aa_substitutions) & !grepl(mut2, aa_substitutions)) |
                                 (grepl(mut2, aa_substitutions) & !grepl(mut1, aa_substitutions))))

#   aa_list <- c("G", "A", "V", "L", "I",
#                "T", "S", "M", "C", "P",
#                "F", "Y", "W", "H", "K",
#                "R", "D", "E", "N", "Q",
#                "stop")
#   mut1 <- row$mutation1
#   mut1_alt <- ifelse(grepl("stop", mut1), 
#                      gsub("stop", "", mut1),
#                      substr(mut1, 1, nchar(mut1) - 1))
#   mut1_alt <- paste0(str_glue("{mut1_alt}{aa_list}"), collapse = "|")
#   
#   mut2 <- row$mutation2
#   mut2_alt <- ifelse(grepl("stop", mut2), 
#                  gsub("stop", "", mut2),
#                  substr(mut2, 1, nchar(mut2) - 1))
#   mut2_alt <- paste0(str_glue("{mut2_alt}{aa_list}"), collapse = "|")
#   
#   res <- gisaid_filt %>%
#     summarise(sum_linked = sum(grepl(mut1, aa_substitutions) & grepl(mut2, aa_substitutions)),
#               sum_unlinked = sum((grepl(mut1, aa_substitutions) & !grepl(mut2_alt, aa_substitutions)) |
#                                    (grepl(mut2, aa_substitutions) & !grepl(mut1_alt, aa_substitutions))))
# }

temp <- row %>%
  mutate(sum_linked = res$sum_linked,
         sum_unlinked = res$sum_unlinked)

fwrite(temp, str_glue("results/allele_frequency_out/linkage_temp/{i}.csv"))

       
  
