rm(list = ls())
setwd("c:/git_repos/wuhu_rooting/")
require(tidyverse)
require(data.table)
require(foreach)
require(randomcoloR)
require(ggrepel)

n <- 1000

alts <- c("A", "T", "G")
ref <- "C"

perm_morsels <- foreach(iter = seq(1000)) %do% {
  morsels <- foreach(i = seq(100), .combine = "c") %do% {
    prob <- runif(1)  
    if(prob >= 0.9) {
        allele <- sample(alts, 1, replace = F)
    } else {
        allele <- ref
    }
    
    return(allele)
  }
  
  tibble(allele = morsels) %>%
    group_by(allele) %>%
    summarise(allele_freq = n() / 1000)
}

bind_rows(perm_morsels) %>%
  filter(allele != ref) %>%
  ggplot(aes(x = allele_freq, color = allele)) +
  geom_density() 
  
