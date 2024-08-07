rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)
require(igraph)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

# Model allele frequencies versus observed 'fitness'
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
monthly_agg <- monthly_df %>%
  group_by(mutation_name, collection_month) %>%
  summarise(monthly_n = sum(n_present),
            monthly_total = unique(n_total)) %>%
  mutate(monthly_prop = monthly_n / monthly_total)

parsed_agg <- monthly_agg %>%
  group_by(mutation_name) %>%
  summarise(global_n = sum(monthly_n),
            max_prop = max(monthly_prop)) %>%
  filter(mutation_name != "")

parsed_filt <- parsed_agg %>% 
  filter(max_prop > 0.1)

obs_df <- fread("results/linkage_out/observed_linkage.csv") %>%
  filter(corr > 0.9) %>%
  filter(mutation1 %in% parsed_filt$mutation_name, 
         mutation2 %in% parsed_filt$mutation_name)

g <- graph_from_data_frame(obs_df %>% dplyr::select(mutation1, mutation2, corr), directed = F)

# g_filt <- delete.edges(g, which(E(g)$corr < 0.9))

l <- layout_with_fr(g)

comm <- cluster_infomap(
  g,e.weights = E(g)$corr
)

plot(comm, g, 
     layout = l,
     vertex.size = 10,
     vertex.label.cex = 0.5)

tibble(mutation1 = names(V(g)),
       comm = comm$membership) %>%
  left_join(obs_df) %>%
  fwrite("results/linkage_out/obs_cor.edgelist.csv")

tibble(mutation1 = names(V(g)),
       comm = comm$membership) %>%
  left_join(obs_df) %>%
  distinct(mutation1) %>%
  nrow()

unique(V(g)$comm)
