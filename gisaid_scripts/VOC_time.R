rm(list = ls())
setwd("c:/git_repos/early_SC2_trajectory/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
require(ghibli)
require(ggpubr)
require(randomcoloR)

gisaid <- fread("../wuhu_rooting/data/metadata/gisaid/gisaid_metatadata.080724.filt.tsv")

monthly_count <- gisaid %>%
  group_by(collection_month) %>%
  summarise(n_total = n())

plot_df <- gisaid %>%
  group_by(collection_month, variant) %>%
  summarise(n = n()) %>%
  mutate(variant_annot = case_when(variant == "" ~ "No designation",
                                   grepl("Alpha", variant) ~ "Alpha",
                                   grepl("Beta", variant) ~ "Beta",
                                   grepl("Gamma", variant) ~ "Gamma",
                                   grepl("Delta", variant) ~ "Delta",
                                   grepl("Omicron", variant) ~ "Omicron")) %>%
  mutate(variant_annot = ifelse(grepl("BA.2.86|JN.1", variant), "BA.2.86/JN.1", variant_annot)) %>%
  mutate(variant_annot = ifelse(grepl("XBB|EG.5", variant), "XBB/EG.5", variant_annot)) %>%
  group_by(variant_annot, collection_month) %>%
  summarise(n = sum(n)) %>%
  left_join(monthly_count) %>%
  mutate(prop = n / n_total) %>%
  mutate(collection_month = as.Date(paste0(collection_month, "-01"))) %>%
  filter(!is.na(variant_annot))


plot_df %>%
  filter(is.na(variant_annot))

# pal <- c(ghibli_palette("KikiMedium")[3:5],
#          ghibli_palette("MarnieMedium1")[7],
#          ghibli_palette("PonyoMedium")[2:3],
#          ghibli_palette("MarnieMedium2")[4],
#          ghibli_palette("MononokeMedium")[3])
# [1] "#5A6F80FF" "#B50A2AFF" "#0E84B4FF" "#E9D097FF" "#9E8356FF" "#278B9AFF" "#44A57CFF"
# [8] "#AE93BEFF"

pal <- c("#5A6F80FF", "#0E84B4FF", "#B50A2AFF","#E9D097FF", "#9E8356FF",
         "#278B9AFF", "#44A57CFF", "#AE93BEFF")

voc_plt <- plot_df %>%
  mutate(variant_annot = factor(variant_annot, c("No designation", "Alpha", "Beta", 
                                                 "Gamma", "Delta", "Omicron",
                                                 "XBB/EG.5", "BA.2.86/JN.1"))) %>%
  ggplot(aes(x = collection_month, y = prop, color = variant_annot)) +
  geom_ribbon(aes(ymin = min(prop),
                  ymax = prop,
                  fill = variant_annot),
              alpha = 0.8) +
  geom_line() +
  geom_point() +
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y",
               limits = c(as.Date("2019-12-01"), as.Date("2024-07-01"))) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal, guide = "none") +
  geom_vline(xintercept = c(as.Date("2021-03-01"), as.Date("2021-06-01"),
                            as.Date("2022-01-01"), as.Date("2023-02-01"),
                            as.Date("2023-12-01")),
             lty = "dashed", color = "black") +
  theme_bw() +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  labs(x = "Collection month", y = "Prop. genomes", fill = "WHO designation")


# gisaid %>% 
#   distinct(variant) %>%
#   mutate(variant_annot = case_when(variant == "" ~ "No designation",
#                                    grepl("Alpha", variant) ~ "Alpha",
#                                    grepl("Beta", variant) ~ "Beta",
#                                    grepl("Gamma", variant) ~ "Gamma",
#                                    grepl("Delta", variant) ~ "Delta",
#                                    grepl("VUM|VOI", variant) ~ "VUM/VOI",
#                                    grepl("Omicron", variant) ~ "Omicron")) %>% View()

monthly_df <- fread("results/allele_frequency_out/all_monthly_frequencies.csv")
parsed_agg <- fread("results/mutation_out/monthly_freq_aggregate.csv")
  
fixed <- parsed_agg %>%
  filter(max_prop > 0.9)

drift <- parsed_agg %>%
  filter(max_prop < 0.1)

major <- parsed_agg %>%
  filter(max_prop > 0.1 & max_prop < 0.9)

interest <- parsed_agg %>%
  filter(mutation_name %in% c("Spike_D614G", "NSP12_P323L", 
                              "N_R203K", "N_G204R", 
                              "NSP4_L438F", "Spike_G142D",
                              "Spike_P681R")) 


mut_pal <- distinctColorPalette(n_distinct(interest$mutation_name))

mut_plt <- monthly_df %>%
  mutate(collection_month = as.Date(paste0(collection_month, "-01"))) %>%
  filter(mutation_name %in% interest$mutation_name) %>%
  ggplot(aes(x = collection_month, y = prop, color = mutation_name, group = mutation_name)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_color_manual(values = mut_pal) +
  geom_vline(xintercept = c(as.Date("2021-03-01"), as.Date("2021-06-01"),
                            as.Date("2022-01-01"), as.Date("2023-02-01"),
                            as.Date("2023-12-01")),
             lty = "dashed", color = "black") +
  theme(text = element_text(family = "sans"),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  scale_x_date(date_breaks = "2 months", 
               date_labels = "%b-%y",
               limits = c(as.Date("2019-12-01"), as.Date("2024-07-01"))) +
  labs(x = "Collection month", y = "Prop. genomes", color = "Mutation")

ggarrange(voc_plt, mut_plt, nrow = 2, align = "hv")

ggsave("results/mutation_out/voc_plot.pdf", dpi = 600, width = 12, height = 4)

