rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.csv")

file_dir <- "results/pipeline_out/stats_out/"
file_list <- list.files(file_dir, full.names = T)
file_list <- file_list[grepl("coverage", file_list)]
length(file_list)

morsels <- foreach(file_name = file_list) %do% {
  id <- gsub(file_dir, "", file_name)
  bs <- gsub(".coverage.txt", "", id)
  
  fread(file_name) %>%
    mutate(biosample = bs) %>%
    mutate(across(everything(), as.character))
}

map_df <- bind_rows(morsels) %>%
  select(biosample, reads_mapped = numreads, coverage_breadth = coverage, mean_depth = meandepth,
         mean_base_quality = meanbaseq, mean_map_quality = meanmapq)

merged <- meta %>%
  left_join(map_df) %>%
  mutate(qc = case_when(is.na(reads_mapped) ~ "sra not processed",
                        coverage_breadth >= 90 & mean_depth >= 20 ~ "high quality",
                        !(coverage_breadth >= 90 & mean_depth >= 20) ~ "low quality"))

# Extract bioprojects that have at least one high quality genome
hq <- merged %>%
  filter(qc == "high quality") 

lq <- merged %>%
  filter(qc == "low quality") %>%
  filter(bioproject_accessions %in% hq$bioproject_accessions)
  
missing <- merged %>%
  filter(qc == "sra not processed") %>%
  filter(bioproject_accessions %in% hq$bioproject_accessions)
 
final <- bind_rows(hq, missing, lq) %>%
  arrange(bioproject_accessions)

final %>%
  ggplot(aes(x = as.Date(collection_date, "%Y-%m-%d"))) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Collection date", y = "No. biosamples")

# Randomly assign bioprojects
set.seed(4)

assignments <- tibble(bioproject_accessions = unique(final$bioproject_accessions),
       curator1 = sample(c(rep("CT", 8),
                  rep("MEZ", 6),
                  rep("FB", 6),
                  rep("AC", 6),
                  rep("SP", 6),
                  rep("SZ", 6),
                  rep("AD", 6),
                  rep("RN", 6),
                  rep("JD", 6),
                  rep("JH", 6),
                  rep("JB", 6)), replace = F))
final %>%
  left_join(assignments) %>%
  group_by(curator1) %>%
  summarise(n_biosamples = n(),
            n_bioprojects = n_distinct(bioproject_accessions))

# Split files
final <- final %>%
  left_join(assignments)

foreach(curator = unique(final$curator1)) %do% {
  final %>%
    filter(curator1 == curator) %>%
    mutate(exclude = NA,
           exclude_reason = NA,
           study_doi = NA,
           hookup_method = NA,
           genome_name_source = NA,
           genome_name = NA,
           gisaid_accession = NA,
           genbank_accession = NA,
           other_accession = NA, .after = 1) %>% 
    mutate(date_verified = NA, .after = 16) %>%
    mutate(loc_verified = NA, .after = 18) %>%
    fwrite(str_glue("data/metadata/sra_metadata/curation_lists/{curator}.csv"))
}
final %>%
  left_join(assignments) %>%
  fwrite("data/metadata/sra_metadata/sra_metadata_curation.180624.csv")

  