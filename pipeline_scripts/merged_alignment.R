rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(Biostrings)

file_list <- list.files("results/pipeline_out/consensus_out/", full.names = T)

ref <- readDNAStringSet("data/genomes/MN908947.3.fna")

fna <- foreach(file_name = file_list, .combine = "c") %do% {
  readDNAStringSet(file_name)  
} 

fna <- c(ref, fna)

writeXStringSet(fna, "data/alignments/reassembled/reassembled.aln")