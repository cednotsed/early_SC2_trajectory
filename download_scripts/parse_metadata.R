rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

gisaid <- fread("data/metadata/gisaid_150324.1Mar20.tsv") %>%
  inner_join(fread("data/metadata/gisaid_150324.1Mar20.more.tsv")) %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  filter(!grepl("Manis|Rhinolo", host)) %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
  separate(location, c(NA, "country"), " / ", remove = F) %>%
  dplyr::rename(accession = accession_id, 
         isolation_source = specimen,
         pango_lineage = lineage,
         release_date = submission_date)
  
ncbi <- fread("data/metadata/ncbi_virus_1Apr20.gt25000.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x)))
  # dplyr::rename(pango_lineage = pangolin)

plot_df <- ncbi %>%
  filter(submitters == "Lv,J.-X., Liu,X., Pei,Y.-Y., Song,Z.-G., Chen,X., Hu,S.-J., Chen,Y.-M., Zhang,Y.-Z.")

plot_df %>%
  ggplot(aes(x = as.Date(collection_date))) +
  geom_histogram() +
  labs(x = "Collection date", y =  "Genomes", title = str_glue("Lv et al. (2024) genomes (n={nrow(plot_df)})"))

# Pango lineages
pango_df <- fread("data/metadata/all_pango_assignments.csv") %>%
  dplyr::rename(accession = taxon) %>%
  separate(qc_notes, c(NA, "prop_ambig"), "\\:") %>%
  select(accession, pango_lineage = lineage, prop_ambig, 
         pango_conflict = conflict, scorpio_call, scorpio_conflict)

pango_df %>%
  filter(!grepl("Omicron|Delta|Alpha|Beta|Gamma|Iota|Theta", scorpio_call)) %>% 
  distinct(scorpio_call)

merged <- bind_rows(gisaid, ncbi) %>%
  separate(collection_date, c("Y", "M", "D"), "-", remove = F) %>%
  select(accession, collection_date, release_date, 
         location, country, host, 
         biosample, bioproject, sra_accession,
         clade) %>%
  left_join(pango_df)
  # Quality filters %>%
  # filter(prop_ambig < 0.05)

merged %>% 
  fwrite("data/metadata/all_genomes.csv", eol = "\n")
# A, B, C
merged %>% 
  filter(grepl("^A\\.", pango_lineage) | pango_lineage == "B"| grepl("^C\\.", pango_lineage)) %>%
  fwrite("data/metadata/parsed_metadata.ABC.csv",
         eol = "\n")

merged %>% 
  filter(country == "China") %>%
  fwrite("data/metadata/parsed_metadata.china_only.csv",
         eol = "\n")

merged %>%
  group_by(country) %>%
  summarise(n = n()) %>% View()
  filter(is.na(M)) %>% 
  group_by(pango_lineage) %>%
  summarise(n = n()) %>% View()

# Get initial complete dataset
# merged %>%
#   filter(pango_lineage == "B.1")
# merged %>%
#   filter(grepl("^A\\.", pango_lineage) | pango_lineage == "B") %>% 
#   filter(is.na(M))

merged %>% 
  filter(!is.na(M) | grepl("^A\\.", pango_lineage) | pango_lineage == "B") %>%
  select(accession, collection_date, release_date, 
         location, country, host, 
         biosample, bioproject, sra_accession,
         pango_lineage, clade) %>%
  fwrite("data/metadata/parsed_metadata.date_complete.csv",
         eol = "\n")

merged %>% 
  filter(!is.na(D)| country == "China") %>% 
  select(accession, collection_date, release_date, 
         location, country, host, 
         biosample, bioproject, sra_accession,
         pango_lineage, clade) %>% 
  fwrite("data/metadata/parsed_metadata.date_complete.csv",
         eol = "\n")
  
merged %>%
  filter(!is.na(D)) %>%
  filter(!is.na(sra_accession)) %>% View()
