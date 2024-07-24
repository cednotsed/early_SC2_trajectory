rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)

gisaid <- fread("data/metadata/gisaid_metatadata.220324.tsv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  rename_all(~gsub("\\?", "", .x)) %>%
  select(accession_id, virus_name, collection_date, 
         location, sequence_length, host, 
         clade, pango_lineage, variant, 
         release_date = submission_date, is_complete, is_high_coverage, 
         is_low_coverage) %>%
  filter(!grepl("Manis|Rhinolo", host)) %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x))) %>%
  separate(location, c(NA, "country"), " / ", remove = F)

gisaid_filt <- gisaid %>%
  filter(collection_date <= "2020-04-01") %>%
  filter(collection_date != "2020") %>%
  filter(is_complete) %>%
  filter(is.na(is_low_coverage)) %>%
  filter(is_high_coverage) %>%
  filter(is.na(variant))

fwrite(gisaid_filt %>% select(virus_name), 
       "data/metadata/gisaid_metatadata.220324.1Apr20.accessions_only.txt",
       eol = "\n",
       col.names = F)

gisaid_filt %>%
  ggplot(aes(x = as.Date(collection_date))) +
  geom_histogram() +
  labs(x = "Collection date", y =  "Genomes", title = str_glue("GISAID genomes up to 1 Apr 2020 (n={nrow(gisaid_filt)})"))
  
group_by(pango_lineage, country) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>% View()


ncbi <- fread("data/metadata/ncbi_virus_1Apr20.gt25000.csv") %>%
  rename_all(~tolower(gsub(" ", "_", .x))) %>%
  mutate(across(everything(), ~ifelse(.x == "", NA, .x)))


ncbi_filt <- ncbi %>%
  filter(length > 29000)

ncbi_filt %>% 
  filter(submitters == "Lv,J.-X., Liu,X., Pei,Y.-Y., Song,Z.-G., Chen,X., Hu,S.-J., Chen,Y.-M., Zhang,Y.-Z.") %>% 
  View()
ncbi_filt %>%
  filter(as.Date(collection_date, "%Y-%m-%d") <= as.Date("2020-03-01")) %>%
  mutate(sra_associated = !is.na(sra_accession)) %>%
  ggplot(aes(x = as.Date(collection_date, "%Y-%m-%d"), fill = sra_associated)) +
  geom_histogram(bins = 100) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", limits = as.Date(c('2019-12-01','2020-04-31'))) +
  labs(x = "Collection date", y = "Genomes", fill = "Has raw reads",
       title = str_glue("NCBI genomes up to 1 Mar 2020")) +
  facet_grid(rows = "sra_associated")

sra_plot <- ncbi_filt %>%
  filter(as.Date(collection_date, "%Y-%m-%d") <= as.Date("2020-03-01")) %>%
  mutate(sra_associated = !is.na(sra_accession)) %>%
  group_by(country) %>%
  summarise(n_sra = sum(sra_associated),
            total = n()) %>%
  mutate(prop = n_sra / total) %>%
  arrange(desc(prop))

sra_plot %>%
  mutate(country = factor(country, unique(sra_plot$country))) %>%
  ggplot(aes(x = country, y = prop, fill = country)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  geom_text(aes(label = n_sra), vjust = 0) +
  labs(x = "Country", y = "Prop. with SRA data")
 
ncbi_filt %>%
  filter(as.Date(collection_date, "%Y-%m-%d") <= as.Date("2020-03-01")) %>%
  mutate(sra_associated = !is.na(sra_accession)) %>%
  filter(sra_associated) %>% View()
ncbi_filt %>% 
  filter(as.Date(collection_date, "%Y-%m-%d") <= as.Date("2020-03-01"))
  filter(!(is.na(biosample) & is.na(bioproject) & is.na(sra_accession))) %>%
  ggplot(aes(x = as.Date(collection_date, "%Y-%m-%d"))) +
  geom_histogram(bins = 100) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", limits = as.Date(c('2019-12-01','2020-04-31'))) +
  labs(x = "Collection date", y = "Genomes", fill = "Lineage",
       title = str_glue("NCBI genomes up to 1 Apr 2020 (n={nrow(ncbi_filt)})")) 
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
