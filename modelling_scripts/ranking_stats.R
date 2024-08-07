rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(Biostrings)
require(foreach)
require(ggpubr)
require(doParallel)
require(relaimpo)

hookup <- fread("data/metadata/SARS-CoV-2_hookup_table_V3.parsed.csv")

df <- fread("results/allele_frequency_out/all_sites.csv") %>%
  left_join(hookup %>% distinct(mutation_name, ref_AA, var_AA))

parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")
monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")

parsed_filt <- parsed_agg %>%
  filter(max_prop > 0.1,
         max_prop <= 0.9)

interest <- df %>%
  arrange(desc(max_freq), desc(n)) %>%
  head(30) %>%
  filter(mutation_name %in% parsed_filt$mutation_name)
  nrow()

  pal <- randomcoloR::distinctColorPalette(n_distinct(interest$mutation_name))
  
  monthly_df %>%
    filter(mutation_name %in% interest$mutation_name) %>%
    mutate(collection_month = as.Date(paste0(collection_month, "-01"))) %>%
    ggplot(aes(x = collection_month, y = prop, color = mutation_name)) +
    geom_point() +
    geom_line() +
    scale_x_date(date_breaks = "2 months", date_labels = "%b-%y") +
    scale_color_manual(values = pal) +
    geom_vline(xintercept = as.Date("2020-03-01"), lty = "dashed") +
    facet_wrap(~mutation_name, ncol = 3) +
    geom_hline(yintercept = 0.1, lty = "dashed") +
    geom_hline(yintercept = 0.9, lty = "dashed") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
          text = element_text(family = "sans")) +
    labs(x = "Collection month", y = "Monthly frequency")

ggsave("results/modelling_out/detected_mutation_frequencies.pdf", dpi = 600, height = 6, width = 10)
  monthly_df %>%
    filter(mutation_name == "Spike_H655Y") %>% View()
  