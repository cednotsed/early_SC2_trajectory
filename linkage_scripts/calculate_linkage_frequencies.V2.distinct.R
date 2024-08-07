rm(list = ls())
setwd("/mnt/c/git_repos/early_SC2_trajectory")
# setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)

args <- commandArgs(trailingOnly = T)
index <- as.numeric(args[1])
# index <- 10

file_dir <- "data/metadata/gisaid/monthly_mutations.distinct/"
file_list <- list.files(file_dir, full.names = T)
file_name <- file_list[index]
# file_name <- "data/metadata/gisaid/monthly_mutations/aa_subs.2022-12.080724.tsv"

# Master df
cor_df <- fread("results/linkage_out/linkage.n5.rsquared25.csv")

# Parse file name
id <- gsub(file_dir, "", file_name)
col_month <- str_split(id, "\\.")[[1]][2]

print(col_month)

set.seed(66)

date_chunk <- fread(file_name, sep = "\t") %>%
  mutate(aa_substitutions = gsub("\\(|\\)", "", aa_substitutions))

n_total <- nrow(date_chunk)

if(n_total > 200000) {
  date_chunk <- date_chunk %>%
    sample_n(200000, replace = F)
}

# Get master list
aa_list <- paste0(date_chunk$aa_substitutions, collapse = ",")
aa_list <- unique(str_split(aa_list, ",")[[1]])

# aa_list <- unique(c(cor_df$mutation1, cor_df$mutation2))

# Create template matrix
mat <- matrix(0, nrow(date_chunk), length(aa_list))
colnames(mat) <- aa_list

# Get presence absence
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

# Filter pairs
cor_filt <- cor_df %>%
  filter(mutation1 %in% aa_list,
         mutation2 %in% aa_list)

# Calculate linked and unlinked counts
link_morsels <- foreach(j = seq(nrow(cor_filt))) %do% {
  if(j %% 10 == 0){
    print(str_glue("{j}/{nrow(cor_filt)}"))
  }
  row <- cor_filt[j, ]
  
  mat_temp <- mat[, unique(c(row$mutation1, row$mutation2))]
  row_sums <- apply(mat_temp, 1, sum)
  n_linked <- sum(row_sums == 2)
  
  # Positive linkage
  if(row$corr > 0) {
    n_unlinked <- sum(row_sums == 1)
  } else {
    # Negative linkage
    possible_variants <- c("G", "A", "V", "L", "I",
                           "T", "S", "M", "C", "P",
                           "F", "Y", "W", "H", "K",
                           "R", "D", "E", "N", "Q",
                           "stop")
    mut1 <- row$mutation1
    mut2 <- row$mutation2
      
    # mut2 and WT 
    mut1_alt <- ifelse(grepl("stop", mut1),
                       gsub("stop", "", mut1),
                       substr(mut1, 1, nchar(mut1) - 1))
    
    mut1_alt <- str_glue("{mut1_alt}{possible_variants}")
    
    alt1 <- colnames(mat)[colnames(mat) %in% mut1_alt]
    
    if(length(alt1) > 0) {
      mat1 <- mat[ , alt1]
      
      if(length(alt1) > 1) {
        mat1 <- rowSums(mat1)
        mat1[mat1 > 1] <- 1
      }
      
      merged_mat1 <- cbind(mat_temp[, 1], mat1)
    } else {
      print("no alt found")
      # If no alts are found, set to zero
      merged_mat1 <- cbind(mat_temp[, 1], rep(0, nrow(mat_temp)))
    }
    
    # mut2 and WT
    mut2_alt <- ifelse(grepl("stop", mut2),
                   gsub("stop", "", mut2),
                   substr(mut2, 1, nchar(mut2) - 1))
    mut2_alt <- str_glue("{mut2_alt}{possible_variants}")
    
    alt2 <- colnames(mat)[colnames(mat) %in% mut2_alt]
    
    if(length(alt2) > 0) {
      mat2 <- mat[ , alt2]
      
      if(length(alt2) > 1) {
        mat2 <- rowSums(mat2)
        mat2[mat2 > 1] <- 1
      }
        
      merged_mat2 <- cbind(mat_temp[, 2], mat2)
    } else {
      # If no alts are found, set to zero
      merged_mat2 <- cbind(mat_temp[, 2], rep(0, nrow(mat_temp)))
    }
    
    union_mat <- rbind(merged_mat1, merged_mat2)
    unlinked_sums <- rowSums(union_mat)
    n_unlinked <- sum(unlinked_sums == 1)
  }
  
  row %>%
    mutate(n_linked = n_linked,
           n_unlinked = n_unlinked)
}

temp <- bind_rows(link_morsels) %>%
  mutate(collection_month = col_month) %>%
  mutate(ratio = n_linked / n_unlinked)

fwrite(temp, str_glue("results/linkage_out/monthly_linkage.distinct/linkage_freq.{col_month}.csv"))

