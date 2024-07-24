rm(list = ls())
setwd("C:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)

meta <- fread("data/metadata/all_metadata.040424.tsv")
pango_meta <- fread("data/metadata/all_metadata.040424.formatted.pango_lineage.csv")

pango_filt <- pango_meta %>% 
  # Retrieve lineage A.*, C.* and B
  filter(grepl("^A\\.|^C\\.", lineage)| lineage == "A"| lineage == "B"| lineage == "C") %>%
  # Remove VoCs and VOIs
  filter(scorpio_call == "") %>%
  select(accession = taxon, pango_lineage = lineage, conflict)

# Filter metadata
meta_filt <- meta %>%
  inner_join(pango_filt)

# Remove sequences with >5% Ns
fna <- readDNAStringSet("data/genomes/all_genomes.040424.formatted.fna")[meta_filt$accession]

parse_sequence <- function(x) {
  seq_mat <- data.matrix(x)
  return(sum(seq_mat == c("N")) / ncol(seq_mat))
}

prop_gaps <- foreach(i = seq(length(fna)), .combine = "bind_rows") %do% {
  tibble(accession = names(fna)[i], prop_gaps = parse_sequence(fna[i]))
}

to_keep <- prop_gaps %>% 
  filter(prop_gaps <= 0.05)

fna_filt <- fna[to_keep$accession]

no_ref <- fna_filt[names(fna_filt) != "MN908947.3"]
ref <- fna_filt[names(fna_filt) == "MN908947.3"]

meta_filt %>%
  filter(accession %in% to_keep$accession) %>%
  fwrite("data/metadata/all_genomes.040424.QCed.ABC.csv")

writeXStringSet(no_ref, "data/genomes/all_genomes.040424.formatted.QCed.ABC.no_ref.fna")
# writeXStringSet(ref, "data/genomes/MN908947.3.fna")
