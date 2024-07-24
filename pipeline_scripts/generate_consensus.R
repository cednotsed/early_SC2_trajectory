rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

non_nanopore <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv") %>%
  filter(platforms != "OXFORD_NANOPORE")

nanopore <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv") %>%
  filter(platforms == "OXFORD_NANOPORE")

ref <- readDNAStringSet("data/genomes/MN908947.3.fna")
ref_vec <- rep("N", width(ref[1]))

file_dir <- "results/pipeline_out/parsed_vcfs/"

file_list <- list.files(file_dir, full.names = T)

foreach(file_name = file_list) %do% {
  # file_name <- file_list[319]
  # file_name = file_list[grepl("SAMN13922059", file_list)]
  id <- gsub(file_dir, "", file_name)
  id <- gsub(".bcftools.parsed.tsv", "", id)
  
  # Coverage file
  cov_df <- fread(str_glue("results/pipeline_out/stats_out/{id}.coverage.txt")) %>%
    mutate(coverage = as.numeric(coverage))
  
  if(cov_df$coverage > 90 & !is.na(cov_df$coverage)) {
    temp <- fread(file_name) %>%
      mutate(across(everything(), as.character)) %>%
      select(-qual) %>%
      separate(ad, c("n_ref", "n_alt1", "n_alt2", "n_alt3"), ",") %>%
      separate(alt, c("alt1", "alt2", "alt3"), ",") %>%
      mutate(across(c(pos, total_depth, mq, n_ref, n_alt1, n_alt2, n_alt3), as.numeric)) %>%
      # Skip indels
      filter(indel_flag == F) %>%
      filter(mq >= 30)
    
    # Non-nanopore filters
    if(id %in% non_nanopore$biosample) {
      temp_filt <- temp %>%
        mutate(consensus = case_when(n_ref > n_alt1 & n_ref >= 10 ~ ref,
                                     n_alt1 > n_ref & n_alt1 >= 10 ~ alt1,
                                     n_alt2 > n_alt1|n_alt3 > n_alt1 ~ "ERROR",
                                     TRUE ~ NA)) %>%
        filter(!is.na(consensus))
    } else if (id %in% nanopore$biosample){
      temp_filt <- temp %>%
        mutate(consensus = case_when(n_ref > n_alt1 & n_ref >= 30 ~ ref,
                                     n_alt1 > n_ref & n_alt1 >= 30 ~ alt1,
                                     n_alt2 > n_alt1|n_alt3 > n_alt1 ~ "ERROR",
                                     TRUE ~ NA)) %>%
        filter(!is.na(consensus))
    } else {
      print(str_glue("{id} is not found in metadata"))
    }
      
    # Check for weird VCFs
    if(any(temp_filt$consensus %in% c("ERROR", "<*>"))) {
      print(str_glue("ERROR in {id}"))
    }
    
    # Write consensus sequence
    consensus <- ref_vec
      
    for(i in seq(nrow(temp_filt))) {
      row <- temp_filt[i, ]
      consensus[row$pos] <- row$consensus
    }
    
    # Print consensus stats
    perc <- sum(consensus == "N") / length(consensus) * 100
    print(str_glue("%Ns = {perc}"))
    
    # Write consensus
    dna <- DNAStringSet(paste0(consensus, collapse = ""))
    names(dna) <- id
    writeXStringSet(dna, str_glue("results/pipeline_out/consensus_out/{id}.fna"))
  }
}  
  

