rm(list = ls())
setwd("C:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(rentrez)

df <- fread("../wuhu_rooting/data/metadata/sra_metadata/biosample_metadata_210324.full.new.tsv")

# df <- df %>%
#   sample_n(1000000, replace = F)

df1 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", taxonomy_name, ignore.case = T))

df2 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", biosample_title, ignore.case = T))

df3 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", isolate, ignore.case = T))

df4 <- df %>%
  filter(grepl("Severe acute respiratory syndrome coronavirus 2|CoV-2|cov2", strain, ignore.case = T))

df_filt <- bind_rows(df1, df2, df3, df4) %>%
  distinct()

dates_filt <- df_filt %>%
  distinct(collection_date) %>% 
  filter(!grepl("/", collection_date)) %>%
  separate(collection_date, c("Y", "M", "D"), "\\-", remove = F) %>%
  filter(!is.na(Y)) %>%
  filter(!is.na(M)) %>% 
  mutate(Y = as.numeric(Y),
         M = as.numeric(M),
         D = as.numeric(D)) %>%
  filter(Y == 2021 & M %in% c(10)|)

final <- df_filt %>%
  filter(collection_date %in% dates_filt$collection_date) %>%
  # Remove samples with "Plasmid" in metadata
  filter(!grepl("plasmid", strain)) %>%
  filter(host != "Plasmid") %>%
  # Remove lab hosts
  filter(lab_host == "") %>%
  # Remove samples with "16S" in metadata
  filter(!grepl("16s|wastewater|environmental", biosample_title)) %>% 
  filter(!grepl("wastewater", isolation_source)) %>% 
  filter(grepl("anterior|nasal|clinical|patient|fecal|human|naso", isolation_source)) %>%
  # Remove samples with missing accessions
  filter(biosample_accession != "") %>%
  left_join(dates_filt)

final %>% View()
  distinct(collection_date, geo_loc_name) %>% View()

# Save biosample accessions to extract SRA run accessions using sra utils
# final %>%
#   dplyr::select(biosample_accession) %>%
#   fwrite("data/metadata/sra_metadata/filtered_sra_accessions.delta.txt", eol = "\n", col.names = F)

# Read run accessions
run_df <- fread("data/metadata/sra_metadata/runinfo.delta.with_headers.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  select(bioproject, biosample, biosample_accession = sample, 
         experiment, run, assemblyname,
         spots, bases, avglength, 
         platform, librarystrategy, libraryselection, 
         librarysource, librarylayout, download_path)

# Get missing efetch queries
final %>%
  left_join(run_df) %>%
  filter(is.na(librarysource)) %>% 
  select(biosample_accession) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.missing.txt", eol = "\n", col.names = F)
# fwrite("data/metadata/sra_metadata/filtered_sra_accessions.csv")

# Parse metadata
parsed <- final %>%
  left_join(run_df) %>%
  select(bioproject, biosample,
         experiment, run, assemblyname,
         collection_date, geo_loc_name, lat_lon,
         host, isolation_source, platform,
         spots, bases, avglength,
         librarystrategy, libraryselection, librarysource,
         librarylayout, download_path)

# Check for double layouts
morsels <- foreach(bs = unique(parsed$biosample)) %do% {
  temp <- parsed %>%
    filter(biosample == bs) %>%
    group_by(biosample, bioproject) %>%
    summarise(n_platform = n_distinct(platform),
              n_layout = n_distinct(librarylayout)) 
  temp
}

dbl_df <- bind_rows(morsels) %>% 
  filter(n_platform > 1| n_layout > 1)

dbl_df %>% 
  group_by(n_platform, n_layout) %>%
  summarise(n = n_distinct(biosample))

# Retrieve runs to delete for double layout biosamples (i.e. retain paired)
to_remove <- parsed %>%
  filter(biosample %in% dbl_df$biosample) %>%
  filter(librarylayout != "PAIRED") %>%
  distinct(run) 

to_remove %>%
  fwrite("data/metadata/sra_metadata/runs_to_remove.txt",
         eol = "\n",
         col.names = F)

parsed_filt <- parsed %>%
  filter(!(run %in% to_remove$run))

# Visualise number of runs
parsed %>% 
  group_by(biosample) %>% 
  summarise(n_runs = n_distinct(run)) %>% 
  group_by(n_runs) %>% 
  summarise(n_biosamples = n()) %>%
  ggplot(aes(x = n_runs, y = n_biosamples)) +
  geom_bar(stat = "identity")
# Merge run accessions
merge_morsels <- foreach(bs = unique(parsed_filt$biosample)) %do% {
  # bs = "SAMN22352829"
  temp <- parsed_filt %>%
    filter(biosample == bs)
  
  dedup <- temp %>%
    distinct(biosample, collection_date, geo_loc_name, host, isolation_source, librarylayout)
  
  if (nrow(dedup) == 1) {
    bases <- temp %>%
      summarise(sum_bases = sum(bases))
    proj_string <- paste0(unique(temp$bioproject), collapse = ";")
    expt_string <- paste0(unique(temp$experiment), collapse = ";")
    run_string <- paste0(unique(temp$run), collapse = ";") 
    platform_string <- paste0(unique(temp$platform, collapse = ";"))
    
    dedup %>%
      mutate(bioproject_accessions = proj_string,
             run_accessions = run_string,
             experiment_accessions = expt_string,
             platforms = platform_string)
  } else {
    print(str_glue("OH NO {bs}")) 
    return(NULL)
  }
}

merged <- bind_rows(merge_morsels) %>%
  relocate(bioproject_accessions, .before = 1) %>%
  relocate(experiment_accessions, .after = 2) %>%
  relocate(run_accessions, platforms, librarylayout, .after = 3) %>%
  mutate(bioproject_accessions = ifelse(bioproject_accessions == "NA", "missing",
                                        bioproject_accessions),
         experiment_accessions = ifelse(experiment_accessions == "NA", "missing",
                                        experiment_accessions),
         run_accessions = ifelse(run_accessions == "NA", "missing",
                                 run_accessions),
         platforms  = ifelse(platforms  == "NA", "missing",
                             platforms),
         librarylayout = ifelse(is.na(librarylayout), "missing",
                                librarylayout)) 

merged %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.csv")

merged %>%
  filter(librarylayout == "PAIRED") %>%
  select(biosample) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.paired.csv",
         eol = "\n",
         col.names = F)

merged %>%
  filter(librarylayout == "SINGLE" & platforms != "OXFORD_NANOPORE") %>%
  select(biosample) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.single.csv",
         eol = "\n",
         col.names = F)

merged %>%
  filter(platforms == "OXFORD_NANOPORE") %>%
  select(biosample) %>%
  fwrite("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.nanopore.csv",
         eol = "\n",
         col.names = F)

final %>% 
  ggplot(aes(x = as.Date(collection_date, "%Y-%m-%d"))) +
  geom_histogram(bins = 100) +
  labs(x = "Collection date", y =  "Genomes", title = str_glue("Biosamples up to 1 Mar 2020 (n={nrow(final)})"))

# df_filt %>%
#   filter(collection_date == "2019-10-26")
# df_filt %>% 
