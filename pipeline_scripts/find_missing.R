setwd("/flask/scratch/matthewsp/early_SC2_trajectory")
require(tidyverse)
require(data.table)
require(foreach)

file_list <- list.files("results/pipeline_out/stats_out", full.names = T)
file_list <- file_list[grepl("coverage", file_list)]
head(file_list)

done <- foreach(file_name = file_list, .combine = "c") %do% {
    print(file_name)
    temp <- fread(file_name)
    id <- gsub("results/pipeline_out/stats_out/|.coverage.txt", "", file_name)

    if(nrow(temp) > 0) {
        #print(id)
        return(id)
    } else {
        #print("NO")
        return(NULL)
    }
}

# NANOPORE
#meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.nanopore.csv", header = F)
#
#parsed <- meta %>%
#    filter(!(V1 %in% done))
#
#print(nrow(parsed))
#
#parsed %>%
#    fwrite("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.nanopore.missing.csv",
#           eol = "\n",
#           col.names = F)
#
# SINGLE
# meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.single.csv", header = F)
#
# parsed <- meta %>%
#     filter(!(V1 %in% done))
#
# print(nrow(parsed))

 #parsed %>%
 #    fwrite("data/metadata/sra_metadata/filtered_sra_accessions.accessions_only.single.missing.csv",
 #           eol = "\n",
 #           col.names = F)

# PAIRED
meta <- fread("data/metadata/sra_metadata/filtered_sra_accessions.transition.accessions_only.paired.csv", header = F)

parsed <- meta %>%
    filter(!(V1 %in% done))

print(nrow(parsed))

 parsed %>%
     fwrite("data/metadata/sra_metadata/filtered_sra_accessions.transition.accessions_only.paired.missing.csv",
            eol = "\n",
            col.names = F)


